#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <HeterogeneousNucleationHandler.h>
#include <experimental/NEReactionNetwork.h>
#include <XolotlConfig.h>
#include <Options.h>
#include <mpi.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace xolotlCore;

class KokkosContext {
public:
	KokkosContext() {
		::Kokkos::initialize();
	}

	~KokkosContext() {
		::Kokkos::finalize();
	}
};
BOOST_GLOBAL_FIXTURE(KokkosContext);

/**
 * This suite is responsible for testing the HeterogeneousNucleationHandler.
 */
BOOST_AUTO_TEST_SUITE(HeterogeneousNucleationHandler_testSuite)

/**
 * Method checking the initialization and the compute re-solution methods.
 */
BOOST_AUTO_TEST_CASE(checkNucleation) {
	// Create the option to create a network
	xolotlCore::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=1000 0 0 0 0" << std::endl << "xenonDiffusivity=1"
			<< std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char **argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array
	// Initialize MPI for HDF5
	MPI_Init(&argc, &argv);
	opts.readParams(argc, argv);

	// Create a grid
	std::vector<double> grid;
	std::vector<double> temperatures;
	int nGrid = 3;
	for (int l = 0; l < nGrid; l++) {
		grid.push_back((double) l * 0.1);
		temperatures.push_back(1800.0);
	}
	// Specify the surface position
	int surfacePos = 0;

	// Create the network
	using NetworkType = experimental::NEReactionNetwork;
	NetworkType::AmountType maxXe = opts.getMaxImpurity();
	NetworkType network( { maxXe }, grid.size(), opts);
	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);
	// Get its size
	const int dof = network.getDOF();

	// Create the re-solution handler
	HeterogeneousNucleationHandler nucleationHandler;

	// Initialize it
	xolotlCore::experimental::IReactionNetwork::SparseFillMap dfill;
	nucleationHandler.initialize(network, dfill);
	nucleationHandler.setFissionYield(0.5);
	nucleationHandler.updateHeterogeneousNucleationRate(1.0);

	// Check some values in dfill
	BOOST_REQUIRE_EQUAL(dfill[0][0], 0);
	BOOST_REQUIRE_EQUAL(dfill[0][1], 1);
	BOOST_REQUIRE_EQUAL(dfill[1][0], 0);
	BOOST_REQUIRE_EQUAL(dfill[1][1], 1);

	// The arrays of concentration
	double concentration[nGrid * dof];
	double newConcentration[nGrid * dof];

	// Initialize their values
	for (int i = 0; i < nGrid * dof; i++) {
		concentration[i] = 110.0;
		newConcentration[i] = 0.0;
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the fifth grid point
	double *concOffset = conc + 1 * dof;
	double *updatedConcOffset = updatedConc + 1 * dof;

	// Set the temperature to compute the rates
	network.setTemperatures(temperatures);
	network.syncClusterDataOnHost();

	// Compute the modified trap mutation at the grid point
	nucleationHandler.computeHeterogeneousNucleation(network, concOffset,
			updatedConcOffset, 1, 0);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 364817.4, 0.01); // Xe
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], -182408.7, 0.01); // Xe_2

	// Initialize the indices and values to set in the Jacobian
	int indices[2];
	double val[2];
	// Get the pointer on them for the compute heterogeneous nucleation method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the heterogeneous nucleation at the grid point 8
	nucleationHandler.computePartialsForHeterogeneousNucleation(network,
			concOffset, valPointer, indicesPointer, 1, 0);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 0); // Xe_1
	BOOST_REQUIRE_EQUAL(indices[1], 1); // Xe_2

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], 3316.5, 0.01); // Xe_1
	BOOST_REQUIRE_CLOSE(val[1], -1658.26, 0.01); // Xe_2

	// Decrease the concentrations to test the regime where there is not enought Xe_1
	for (int i = 0; i < nGrid * dof; i++) {
		concentration[i] = 50.0;
		newConcentration[i] = 0.0;
	}

	// Compute the modified trap mutation at the grid point
	nucleationHandler.computeHeterogeneousNucleation(network, concOffset,
			updatedConcOffset, 1, 0);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 75348.22, 0.01); // Xe
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], -37674.11, 0.01); // Xe_2

	// Compute the partial derivatives for the heterogeneous nucleation at the grid point 8
	nucleationHandler.computePartialsForHeterogeneousNucleation(network,
			concOffset, valPointer, indicesPointer, 1, 0);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 0); // Xe_1
	BOOST_REQUIRE_EQUAL(indices[1], 1); // Xe_2

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], 1506.96, 0.01); // Xe_1
	BOOST_REQUIRE_CLOSE(val[1], -753.48, 0.01); // Xe_2

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
