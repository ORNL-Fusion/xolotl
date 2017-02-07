#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <DummyTrapMutationHandler.h>
#include <Sigma3TrapMutationHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <DummyAdvectionHandler.h>
#include <YGBAdvectionHandler.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the Sigma3TrapMutationHandler.
 */
BOOST_AUTO_TEST_SUITE(Sigma3TrapMutationHandler_testSuite)

/**
 * Method checking the initialization and the compute modified trap-mutation methods.
 */
BOOST_AUTO_TEST_CASE(checkModifiedTrapMutation) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = loader.load().get();
	// Get all the reactants
	auto allReactants = network->getAll();
	// Get its size
	const int size = network->size();
	const int dof = network->getDOF();
	// Initialize the rate constants
	for (int i = 0; i < size; i++) {
		// This part will set the temperature in each reactant
		// and recompute the diffusion coefficient
		allReactants->at(i)->setTemperature(1000.0);
	}
	for (int i = 0; i < size; i++) {
		// Now that the diffusion coefficients of all the reactants
		// are updated, the reaction and dissociation rates can be
		// recomputed
		auto cluster = (xolotlCore::PSICluster *) allReactants->at(i);
		cluster->computeRateConstants();
	}

	// Suppose we have a grid with 13 grip points and distance of
	// 0.1 nm between grid points
	std::vector<double> grid;
	for (int l = 0; l < 13; l++) {
		grid.push_back((double) l * 0.1);
	}
	// Set the surface position
	std::vector<int> surfacePos = { 0, 0, 0, 0, 0 };

	// Create the modified trap-mutation handler
	DummyTrapMutationHandler trapMutationHandler;

	// Create the advection handlers needed to initialize the trap mutation handler
	std::vector<xolotlCore::IAdvectionHandler *> advectionHandlers;
	advectionHandlers.push_back(new DummyAdvectionHandler());
	auto advecHandler = new YGBAdvectionHandler();
	advecHandler->setLocation(1.0);
	advecHandler->setDimension(2);
	advectionHandlers.push_back(advecHandler);

	// Initialize it
	trapMutationHandler.initialize(network, grid, 5, 0.5);
	trapMutationHandler.initializeIndex2D(surfacePos, network,
			advectionHandlers, grid, 5, 0.5);

	// The arrays of concentration
	double concentration[13 * 5 * dof];
	double newConcentration[13 * 5 * dof];

	// Initialize their values
	for (int i = 0; i < 13 * 5 * dof; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the sixth grid point on the second row
	double *concOffset = conc + (13 * 1 + 5) * dof;
	double *updatedConcOffset = updatedConc + (13 * 1 + 5) * dof;

	// Putting the concentrations in the network so that the rate for
	// desorption is computed correctly
	network->updateConcentrationsFromArray(concOffset);

	// Compute the modified trap mutation at the sixth grid point
	trapMutationHandler.computeTrapMutation(network, concOffset,
			updatedConcOffset, 5, 1);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 5.35980e+23, 0.01); // Create I
	BOOST_REQUIRE_CLOSE(updatedConcOffset[9], -1.33984e+23, 0.01); // He4
	BOOST_REQUIRE_CLOSE(updatedConcOffset[18], 1.33984e+23, 0.01); // Create He4V

	// Get the offset for the ninth grid point on the fourth row
	concOffset = conc + (13 * 3 + 8) * dof;
	updatedConcOffset = updatedConc + (13 * 3 + 8) * dof;

	// Putting the concentrations in the network so that the rate for
	// desorption is computed correctly
	network->updateConcentrationsFromArray(concOffset);

	// Compute the modified trap mutation at the ninth grid point
	trapMutationHandler.computeTrapMutation(network, concOffset,
			updatedConcOffset, 8, 3);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 3.65299e+24, 0.01); // Create I
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8], 0.0, 0.01); // He3
	BOOST_REQUIRE_CLOSE(updatedConcOffset[17], 0.0, 0.01); // Doesn't create He3V
	BOOST_REQUIRE_CLOSE(updatedConcOffset[12], -9.132756e+23, 0.01); // He7
	BOOST_REQUIRE_CLOSE(updatedConcOffset[21], 9.132756e+23, 0.01); // Create He7V2

	// Initialize the indices and values to set in the Jacobian
	int nHelium = network->getAll(heType).size();
	int indices[3 * nHelium];
	double val[3 * nHelium];
	// Get the pointer on them for the compute modified trap-mutation method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the modified trap-mutation at the grid point 8
	int nMutating = trapMutationHandler.computePartialsForTrapMutation(network,
			valPointer, indicesPointer, 8, 3);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(nMutating, 4);
	BOOST_REQUIRE_EQUAL(indices[0], 9); // He4
	BOOST_REQUIRE_EQUAL(indices[1], 18); // He4V
	BOOST_REQUIRE_EQUAL(indices[2], 0); // I
	BOOST_REQUIRE_EQUAL(indices[3], 10); // He5
	BOOST_REQUIRE_EQUAL(indices[4], 19); // He5V
	BOOST_REQUIRE_EQUAL(indices[5], 0); // I

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], -9.67426e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 9.67426e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 9.67426e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], -9.67426e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], 9.67426e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[5], 9.67426e+13, 0.01);

	// Change the temperature of the network
	network->setTemperature(500.0);

	// Update the bursting rate
	trapMutationHandler.updateTrapMutationRate(network);

	// Compute the partial derivatives for the bursting a the grid point 8
	nMutating = trapMutationHandler.computePartialsForTrapMutation(network,
			valPointer, indicesPointer, 8, 3);

	// Check values
	BOOST_REQUIRE_EQUAL(nMutating, 4);
	BOOST_REQUIRE_CLOSE(val[0], -2.14016e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 2.14016e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 2.14016e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], -2.14016e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], 2.14016e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[5], 2.14016e+13, 0.01);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
