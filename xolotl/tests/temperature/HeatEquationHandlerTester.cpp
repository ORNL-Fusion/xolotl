#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <HeatEquationHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <Options.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the HeatEquationHandler.
 */
BOOST_AUTO_TEST_SUITE(Diffusion1DHandler_testSuite)

/**
 * Method checking the initialization of the off-diagonal and diagonal part of the Jacobian,
 * and the compute temperature methods.
 */
BOOST_AUTO_TEST_CASE(checkHeat) {
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
	for (int l = 0; l < 3; l++) {
		grid.push_back((double) l);
	}

	// Create the heat handler
	HeatEquationHandler heatHandler = HeatEquationHandler(1200.0, 1000.0);
	heatHandler.setHeatCoefficient(xolotlCore::tungstenHeatCoefficient);

	// Check the initial temperatures
	BOOST_REQUIRE_CLOSE(heatHandler.getTemperature( { 0.0, 0.0, 0.0 }, 0.0),
			1200.0, 0.01);
	BOOST_REQUIRE_CLOSE(heatHandler.getTemperature( { 1.0, 0.0, 0.0 }, 0.0),
			1000.0, 0.01);

	// Create ofill
	int matO[dof * dof];
	int *ofill = &matO[0];
	// Create dfill
	int matD[dof * dof];
	int *dfill = &matD[0];

	// Initialize it
	heatHandler.initializeTemperature(*network, ofill, dfill);

	// Check that the temperature "diffusion" is well set
	BOOST_REQUIRE_EQUAL(ofill[99], 1);
	BOOST_REQUIRE_EQUAL(dfill[99], 1);

	// The size parameter in the x direction
	double hx = 1.0;

	// The arrays of concentration
	double concentration[3 * dof];
	double newConcentration[3 * dof];

	// Initialize their values
	for (int i = 0; i < 3 * dof; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	// Supposing the 3 grid points are laid-out as follow:
	// 0 | 1 | 2
	double *concOffset = conc + dof;
	double *updatedConcOffset = updatedConc + dof;

	// Fill the concVector with the pointer to the middle, left, and right grid points
	double **concVector = new double*[3];
	concVector[0] = concOffset; // middle
	concVector[1] = conc; // left
	concVector[2] = conc + 2 * dof; // right

	// Compute the heat equation at this grid point
	heatHandler.computeTemperature(concVector, updatedConcOffset, hx, hx);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[9], 1.367e+16, 0.01);

	// Set the temperature in the handler
	heatHandler.setTemperature(concOffset);
	// Check the updated temperature
	Point<3> pos { 1.0, 0.0, 0.0 };
	BOOST_REQUIRE_CLOSE(heatHandler.getTemperature(pos, 1.0), 361.0, 0.01);

	// Initialize the indices and values to set in the Jacobian
	int indices[1];
	double val[3];
	// Get the pointer on them for the compute diffusion method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the heat equation a the grid point
	heatHandler.computePartialsForTemperature(valPointer, indicesPointer, hx,
			hx);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 9);

	// Check the values
	BOOST_REQUIRE_CLOSE(val[0], -1.367e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 6.835e+13, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 6.835e+13, 0.01);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
