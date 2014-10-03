#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <DiffusionHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the DiffusionHandler.
 */
BOOST_AUTO_TEST_SUITE(DiffusionHandler_testSuite)

/**
 * Method checking the initialization of the off-diagonal part of the Jacobian,
 * and the compute diffusion methods.
 */
BOOST_AUTO_TEST_CASE(checkDiffusion) {
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

	// Create the diffusion handler
	DiffusionHandler diffusionHandler;

	// Create ofill
	PetscInt mat[size*size];
	PetscInt *ofill = &mat[0];

	// Initialize it
	diffusionHandler.initializeOFill(network, ofill);

	// All the clusters diffuse except the last one
	BOOST_REQUIRE_EQUAL(ofill[0], 1);
	BOOST_REQUIRE_EQUAL(ofill[6], 1);
	BOOST_REQUIRE_EQUAL(ofill[12], 1);
	BOOST_REQUIRE_EQUAL(ofill[18], 1);

	// Check the total number of diffusing clusters
	BOOST_REQUIRE_EQUAL(diffusionHandler.getNumberOfDiffusing(), 4);

	// The size parameter
	double sx = 1.0;

	// The arrays of concentration
	PetscReal concentration[3*size];
	PetscReal newConcentration[3*size];

	// Initialize their values
	for (int i = 0; i < 3*size; i++) {
		concentration[i] = (PetscReal) i * i;
		newConcentration[i] = 0.0;
	}

	// Set the temperature to 1000 K to initialize the diffusion coefficients
	auto reactants = network->getAll();
	for (int i = 0; i < size; i++) {
		auto cluster = (PSICluster *) reactants->at(i);
		cluster->setTemperature(1000.0);
	}

	// Get pointers
	PetscReal *conc = &concentration[0];
	PetscReal *updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	PetscReal *concOffset = conc + size ;
	PetscReal *leftConcOffset = conc;
	PetscReal *rightConcOffset = conc + size * 2;
	PetscReal *updatedConcOffset = updatedConc + size;

	// Compute the diffusion at this grid point
	diffusionHandler.computeDiffusion(network, sx, concOffset,
			leftConcOffset, rightConcOffset, updatedConcOffset);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 3917913916907.0, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], 320772236824.0, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[2], 157095680823.0, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], 25265634.0, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4], 0.0, 0.01); // Does not diffuse

	// Initialize the rows, columns, and values to set in the Jacobian
	int nDiff = diffusionHandler.getNumberOfDiffusing();
	PetscInt row[nDiff], col[3*nDiff];
	PetscReal val[3*nDiff];
	// Get the pointer on them for the compute diffusion method
	PetscInt *rowPointer = &row[0];
	PetscInt *colPointer = &col[0];
	PetscReal *valPointer = &val[0];

	// Compute the partial derivatives for the diffusion a the grid point 1
	diffusionHandler.computePartialsForDiffusion(network, sx, valPointer,
			rowPointer, colPointer, 1, 0);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(row[0], 10);
	BOOST_REQUIRE_EQUAL(row[1], 11);
	BOOST_REQUIRE_EQUAL(row[2], 12);
	BOOST_REQUIRE_EQUAL(row[3], 13);

	BOOST_REQUIRE_EQUAL(col[0], 5);
	BOOST_REQUIRE_EQUAL(col[1], 10);
	BOOST_REQUIRE_EQUAL(col[2], 15);
	BOOST_REQUIRE_EQUAL(col[3], 6);
	BOOST_REQUIRE_EQUAL(col[4], 11);
	BOOST_REQUIRE_EQUAL(col[5], 16);
	BOOST_REQUIRE_EQUAL(col[6], 7);
	BOOST_REQUIRE_EQUAL(col[7], 12);
	BOOST_REQUIRE_EQUAL(col[8], 17);
	BOOST_REQUIRE_EQUAL(col[9], 8);
	BOOST_REQUIRE_EQUAL(col[10], 13);
	BOOST_REQUIRE_EQUAL(col[11], 18);

	// Check some values
	BOOST_REQUIRE_CLOSE(val[0], 78358278338.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 6415444736.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[5], 6415444736.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[7], -6283827232.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[10], -1010625.0, 0.01);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
