#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <ZGBAdvectionHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the ZGBAdvectionHandler.
 */
BOOST_AUTO_TEST_SUITE(ZGBAdvectionHandler_testSuite)

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

	// Load the network
	auto network = loader.load().get();
	// Get its size
	const int dof = network->getDOF();

	// Create ofill
	int mat[dof * dof];
	int *ofill = &mat[0];

	// Create the advection handler and initialize it with a sink at
	// 2nm in the Z direction
	ZGBAdvectionHandler advectionHandler;
	advectionHandler.initialize(network, ofill);
	advectionHandler.setLocation(2.0);
	advectionHandler.setDimension(3);

	// Check if grid points are on the sink
	std::vector<double> pos0 = { 0.1, 3.0, 0.0 };
	std::vector<double> pos1 = { 2.0, 2.0, 2.0 };
	BOOST_REQUIRE_EQUAL(advectionHandler.isPointOnSink(pos0), false);
	BOOST_REQUIRE_EQUAL(advectionHandler.isPointOnSink(pos1), true);

	// Check the total number of advecting clusters
	BOOST_REQUIRE_EQUAL(advectionHandler.getNumberOfAdvecting(), 6);

	// Set the size parameters
	double hx = 1.0;
	double hy = 0.5;
	double hz = 2.0;

	// The arrays of concentration
	double concentration[27 * dof];
	double newConcentration[27 * dof];

	// Initialize their values
	for (int i = 0; i < 27 * dof; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Set the temperature to 1000K to initialize the diffusion coefficients
	auto reactants = network->getAll();
	for (int i = 0; i < dof - 1; i++) {
		auto cluster = (PSICluster *) reactants->at(i);
		cluster->setTemperature(1000.0);
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	// Supposing the 27 grid points are laid-out as follow (a cube!):
	// 6 | 7 | 8    15 | 16 | 17    24 | 25 | 26
	// 3 | 4 | 5    12 | 13 | 14    21 | 22 | 23
	// 0 | 1 | 2    9  | 10 | 11    18 | 19 | 20
	//   front         middle           back
	double *concOffset = conc + 13 * dof;
	double *updatedConcOffset = updatedConc + 13 * dof;

	// Fill the concVector with the pointer to the middle, left, right, bottom, top, front, and back grid points
	double **concVector = new double*[7];
	concVector[0] = concOffset; // middle
	concVector[1] = conc + 12 * dof; // left
	concVector[2] = conc + 14 * dof; // right
	concVector[3] = conc + 10 * dof; // bottom
	concVector[4] = conc + 16 * dof; // top
	concVector[5] = conc + 4 * dof; // front
	concVector[6] = conc + 22 * dof; // back

	// Set the grid position
	std::vector<double> gridPosition = { hx, hy, hz };

	// Compute the advection at this grid point
	advectionHandler.computeAdvection(network, gridPosition, concVector,
			updatedConcOffset, hx, hx, 1, hy, 1, hz, 1);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 1.88447e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], 1.74419e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[2], 2.12717e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], 3.68085e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4], 5.09884e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[5], 1.90674e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[6], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[7], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8], 0.0, 0.01); // Does not advect

	// Initialize the rows, columns, and values to set in the Jacobian
	int nAdvec = advectionHandler.getNumberOfAdvecting();
	int indices[nAdvec];
	double val[7 * nAdvec];
	// Get the pointer on them for the compute advection method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the advection a the grid point 1
	advectionHandler.computePartialsForAdvection(network, valPointer,
			indicesPointer, gridPosition, hx, hx, 1, hy, 1, hz, 1);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 0);
	BOOST_REQUIRE_EQUAL(indices[1], 1);
	BOOST_REQUIRE_EQUAL(indices[2], 2);
	BOOST_REQUIRE_EQUAL(indices[3], 3);
	BOOST_REQUIRE_EQUAL(indices[4], 4);
	BOOST_REQUIRE_EQUAL(indices[5], 5);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], 3.76893e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 3.76893e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 3.45234e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 3.45234e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], 4.16701e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[5], 4.16701e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[6], 7.13648e+06, 0.01);

	// Get the stencil
	auto stencil = advectionHandler.getStencilForAdvection(gridPosition);

	// Check the value of the stencil
	BOOST_REQUIRE_EQUAL(stencil[0], 0);
	BOOST_REQUIRE_EQUAL(stencil[1], 0);
	BOOST_REQUIRE_EQUAL(stencil[2], 1); // z

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
