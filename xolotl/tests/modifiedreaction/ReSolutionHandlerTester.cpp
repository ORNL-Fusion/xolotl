#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <ReSolutionHandler.h>
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
 * This suite is responsible for testing the ReSolutionHandler.
 */
BOOST_AUTO_TEST_SUITE(ReSolutionHandler_testSuite)

/**
 * Method checking the initialization and the compute re-solution methods.
 */
BOOST_AUTO_TEST_CASE(checkReSolution) {
	// Create the option to create a network
	xolotlCore::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=10000 0 0 0 0" << std::endl;
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
	// Initialize MPI
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
	ReSolutionHandler reSolutionHandler;

	// Initialize it
	xolotlCore::experimental::IReactionNetwork::SparseFillMap dfill;
	reSolutionHandler.initialize(network, dfill, 0.73);
	reSolutionHandler.setFissionYield(0.25);
	reSolutionHandler.updateReSolutionRate(1.0);

	// Check some values in dfill
	BOOST_REQUIRE_EQUAL(dfill[1][0], 1);
	BOOST_REQUIRE_EQUAL(dfill[3][0], 3);
	BOOST_REQUIRE_EQUAL(dfill[5][0], 5);
	BOOST_REQUIRE_EQUAL(dfill[7][0], 7);
	BOOST_REQUIRE_EQUAL(dfill[9][0], 9);
	BOOST_REQUIRE_EQUAL(dfill[11][0], 11);
	BOOST_REQUIRE_EQUAL(dfill[13][0], 13);

	// The arrays of concentration
	double concentration[nGrid * dof];
	double newConcentration[nGrid * dof];

	// Initialize their values
	for (int i = 0; i < nGrid * dof; i++) {
		concentration[i] = (double) i * i;
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

	// Compute the modified trap mutation at the sixth grid point
	reSolutionHandler.computeReSolution(network, concOffset, updatedConcOffset,
			1, 0);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 6.5039236e+16, 0.01); // Create Xe
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8000], 7.6687556e+8, 0.01); // Xe_7999
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8001], 7.66911126e+8, 0.01); // Xe_8000

	// Initialize the indices and values to set in the Jacobian
	int nXenon = reSolutionHandler.getNumberOfReSoluting();
	int indices[10 * nXenon];
	double val[10 * nXenon];
	// Get the pointer on them for the compute re-solution method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the re-solution at the grid point 8
	int nReSo = reSolutionHandler.computePartialsForReSolution(network,
			valPointer, indicesPointer, 1, 0);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(nReSo, 9999);
	BOOST_REQUIRE_EQUAL(indices[0], 1); // Xe_2
	BOOST_REQUIRE_EQUAL(indices[1], 1); // Xe_2
	BOOST_REQUIRE_EQUAL(indices[2], 0); // Xe_1
	BOOST_REQUIRE_EQUAL(indices[3], 0); // Xe_1
	BOOST_REQUIRE_EQUAL(indices[4], 0); // Xe_1

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], -197711.9, 0.01); // Xe_2
	BOOST_REQUIRE_CLOSE(val[1], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[2], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[3], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[4], 197711.9, 0.01); // Xe_1
	BOOST_REQUIRE_CLOSE(val[5], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[6], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[7], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[8], 197711.9, 0.01); // Xe_1
	BOOST_REQUIRE_CLOSE(val[9], 0.0, 0.01); // no grouping

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * Method checking the use of a minimum size.
 */
BOOST_AUTO_TEST_CASE(checkMinimumSize) {
	// Create the option to create a network
	xolotlCore::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=10000 0 0 0 0" << std::endl;
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
	ReSolutionHandler reSolutionHandler;

	// Initialize it
	xolotlCore::experimental::IReactionNetwork::SparseFillMap dfill;
	reSolutionHandler.initialize(network, dfill, 0.73);
	reSolutionHandler.setFissionYield(0.25);
	reSolutionHandler.updateReSolutionRate(1.0);
	reSolutionHandler.setMinSize(10);

	// The arrays of concentration
	double concentration[nGrid * dof];
	double newConcentration[nGrid * dof];

	// Initialize their values
	for (int i = 0; i < nGrid * dof; i++) {
		concentration[i] = (double) i * i;
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

	// Compute the modified trap mutation at the sixth grid point
	reSolutionHandler.computeReSolution(network, concOffset, updatedConcOffset,
			1, 0);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 6.503923568e+16, 0.01); // Create Xe
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8000], 766875563, 0.01); // Xe_7999
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8001], 766875563, 0.01); // Xe_8000

	// Initialize the indices and values to set in the Jacobian
	int nXenon = reSolutionHandler.getNumberOfReSoluting();
	int indices[10 * nXenon];
	double val[10 * nXenon];
	// Get the pointer on them for the compute re-solution method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the re-solution at the grid point 8
	int nReSo = reSolutionHandler.computePartialsForReSolution(network,
			valPointer, indicesPointer, 1, 0);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], -197711.9, 0.01); // Xe_2
	BOOST_REQUIRE_CLOSE(val[1], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[2], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[3], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[4], 197711.9, 0.01); // Xe_1
	BOOST_REQUIRE_CLOSE(val[5], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[6], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[7], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[8], 197711.9, 0.01); // Xe_1
	BOOST_REQUIRE_CLOSE(val[9], 0.0, 0.01); // no grouping

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * Method checking the use of a different fit.
 */
BOOST_AUTO_TEST_CASE(checkDifferentFit) {
	// Create the option to create a network
	xolotlCore::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=10000 0 0 0 0" << std::endl;
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
	ReSolutionHandler reSolutionHandler;

	// Initialize it
	xolotlCore::experimental::IReactionNetwork::SparseFillMap dfill;
	reSolutionHandler.initialize(network, dfill, 1.0);
	reSolutionHandler.setFissionYield(0.25);
	reSolutionHandler.updateReSolutionRate(1.0);

	// The arrays of concentration
	double concentration[nGrid * dof];
	double newConcentration[nGrid * dof];

	// Initialize their values
	for (int i = 0; i < nGrid * dof; i++) {
		concentration[i] = (double) i * i;
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

	// Compute the modified trap mutation at the sixth grid point
	reSolutionHandler.computeReSolution(network, concOffset, updatedConcOffset,
			1, 0);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 1.2634756e+17, 0.01); // Create Xe
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8000], 1.609505e+9, 0.01); // Xe_7999
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8001], 1.609589e+9, 0.01); // Xe_8000

	// Initialize the indices and values to set in the Jacobian
	int nXenon = reSolutionHandler.getNumberOfReSoluting();
	int indices[10 * nXenon];
	double val[10 * nXenon];
	// Get the pointer on them for the compute re-solution method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the re-solution at the grid point 8
	int nReSo = reSolutionHandler.computePartialsForReSolution(network,
			valPointer, indicesPointer, 1, 0);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], -324044, 0.01); // Xe_2
	BOOST_REQUIRE_CLOSE(val[1], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[2], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[3], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[4], 324044, 0.01); // Xe_1
	BOOST_REQUIRE_CLOSE(val[5], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[6], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[7], 0.0, 0.01); // no grouping
	BOOST_REQUIRE_CLOSE(val[8], 324044, 0.01); // Xe_1
	BOOST_REQUIRE_CLOSE(val[9], 0.0, 0.01); // no grouping

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
