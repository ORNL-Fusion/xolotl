#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <W111AdvectionHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the W111AdvectionHandler.
 */
BOOST_AUTO_TEST_SUITE(W111AdvectionHandler_testSuite)

/**
 * Method checking the initialization and the compute advection methods.
 */
BOOST_AUTO_TEST_CASE(checkAdvection) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = loader.load();
	// Get its size
	const int size = network->getAll()->size();

	// Create the advection handler
	W111AdvectionHandler advectionHandler;

	// Initialize it
	advectionHandler.initialize(network);

	// Check the total number of advecting clusters
	BOOST_REQUIRE_EQUAL(advectionHandler.getNumberOfAdvecting(), 2);

	// The size parameter
	double hx = 1.0;

	// The arrays of concentration
	double concentration[3*size];
	double newConcentration[3*size];

	// Initialize their values
	for (int i = 0; i < 3*size; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Set the temperature to 1000 K to initialize the diffusion coefficients
	auto reactants = network->getAll();
	for (int i = 0; i < size; i++) {
		auto cluster = (PSICluster *) reactants->at(i);
		cluster->setTemperature(1000.0);
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	double *concOffset = conc + size;
	double *rightConcOffset = conc + size * 2;
	double *updatedConcOffset = updatedConc + size;

	// Compute the advection at this grid point
	advectionHandler.computeAdvection(network, hx, 1,
			concOffset, rightConcOffset, updatedConcOffset);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], -23182456628.0, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[2], -28001585008.0, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4], 0.0, 0.01); // Does not advect

	// Initialize the rows, columns, and values to set in the Jacobian
	int nAdvec = advectionHandler.getNumberOfAdvecting();
	int row[nAdvec], col[2*nAdvec];
	double val[2*nAdvec];
	// Get the pointer on them for the compute advection method
	int *rowPointer = &row[0];
	int *colPointer = &col[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the advection a the grid point 1
	advectionHandler.computePartialsForAdvection(network, hx, valPointer,
			rowPointer, colPointer, 1, 0);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(row[0], 11);
	BOOST_REQUIRE_EQUAL(row[1], 12);

	BOOST_REQUIRE_EQUAL(col[0], 6);
	BOOST_REQUIRE_EQUAL(col[1], 11);
	BOOST_REQUIRE_EQUAL(col[2], 7);
	BOOST_REQUIRE_EQUAL(col[3], 12);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], 815207266.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], -815207266.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 700039625.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], -700039625.0, 0.01);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
