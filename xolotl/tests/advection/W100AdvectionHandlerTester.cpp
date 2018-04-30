#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <W100AdvectionHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <Options.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the W100AdvectionHandler.
 */
BOOST_AUTO_TEST_SUITE(W100AdvectionHandler_testSuite)

/**
 * Method checking the initialization and the compute advection methods.
 */
BOOST_AUTO_TEST_CASE(checkAdvection) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);
	// Get its size
	const int dof = network->getDOF();

	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double) l);
	}

	// Create ofill
	int mat[dof * dof];
	int *ofill = &mat[0];

	// Create a collection of advection handlers
	std::vector<IAdvectionHandler *> advectionHandlers;

	// Create the advection handler and initialize it
	W100AdvectionHandler advectionHandler;
	advectionHandler.initialize(*network, ofill);
	advectionHandler.initializeAdvectionGrid(advectionHandlers, grid);

	// Check the total number of advecting clusters
	BOOST_REQUIRE_EQUAL(advectionHandler.getNumberOfAdvecting(), 7);

	// Check the clusters in ofill
	BOOST_REQUIRE_EQUAL(ofill[0], 1);
	BOOST_REQUIRE_EQUAL(ofill[11], 1);
	BOOST_REQUIRE_EQUAL(ofill[22], 1);
	BOOST_REQUIRE_EQUAL(ofill[33], 1);
	BOOST_REQUIRE_EQUAL(ofill[44], 1);
	BOOST_REQUIRE_EQUAL(ofill[55], 1);
	BOOST_REQUIRE_EQUAL(ofill[66], 1);

	// Set the size parameter in the x direction
	double hx = 1.0;

	// Create the arrays of concentration
	double concentration[3 * dof];
	double newConcentration[3 * dof];

	// Initialize their values
	for (int i = 0; i < 3 * dof; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Set the temperature to 1000K to initialize the diffusion coefficients
	network->setTemperature(1000.0);

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	double *concOffset = conc + dof;
	double *updatedConcOffset = updatedConc + dof;

	// Fill the concVector with the pointer to the middle, left, and right grid points
	double **concVector = new double*[3];
	concVector[0] = concOffset; // middle
	concVector[1] = conc; // left
	concVector[2] = conc + 2 * dof; // right

	// Set the grid position
	Point<3> gridPosition { hx, 0.0, 0.0 };

	// Compute the advection at this grid point
	advectionHandler.computeAdvection(*network, gridPosition, concVector,
			updatedConcOffset, hx, hx, 1);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], -3.81919e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], -5.17147e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[2], -3.63429e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], -1.25359e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4], -1.17283e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[5], -5.40893e+09, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[6], -1.27548e+09, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[7], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8], 0.0, 0.01); // Does not advect

	// Initialize the rows, columns, and values to set in the Jacobian
	int nAdvec = advectionHandler.getNumberOfAdvecting();
	int indices[nAdvec];
	double val[2 * nAdvec];
	// Get the pointer on them for the compute advection method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the advection a the grid point 1
	advectionHandler.computePartialsForAdvection(*network, valPointer,
			indicesPointer, gridPosition, hx, hx, 1);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 0);
	BOOST_REQUIRE_EQUAL(indices[1], 1);
	BOOST_REQUIRE_EQUAL(indices[2], 2);
	BOOST_REQUIRE_EQUAL(indices[3], 3);
	BOOST_REQUIRE_EQUAL(indices[4], 4);
	BOOST_REQUIRE_EQUAL(indices[5], 5);
	BOOST_REQUIRE_EQUAL(indices[6], 6);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], -509225360.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 31826585.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], -553468828.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 34591802.0, 0.01);

	// Get the stencil
	auto stencil = advectionHandler.getStencilForAdvection(gridPosition);

	// Check the value of the stencil
	BOOST_REQUIRE_EQUAL(stencil[0], 1); //x
	BOOST_REQUIRE_EQUAL(stencil[1], 0);
	BOOST_REQUIRE_EQUAL(stencil[2], 0);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
