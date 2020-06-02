#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <xolotl/core/flux/W211FitFluxHandler.h>
#include <xolotl/config.h>
#include <xolotl/options/Options.h>

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
 * The test suite is responsible for testing the W211FitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (W211FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkComputeIncidentFlux) {
	// Create the option to create a network
	xolotlCore::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=9 0 0 0 0" << std::endl;
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
	for (int l = 0; l < 7; l++) {
		grid.push_back((double) l * 1.25);
	}
	// Specify the surface position
	int surfacePos = 0;

	// Create the network
	using NetworkType =
	experimental::PSIReactionNetwork<experimental::PSIFullSpeciesList>;
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe = opts.getMaxImpurity();
	NetworkType::AmountType maxD = opts.getMaxD();
	NetworkType::AmountType maxT = opts.getMaxT();
	NetworkType network( { maxHe, maxD, maxT, maxV, maxI }, grid.size(), opts);
	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);
	// Get its size
	const int dof = network.getDOF();

	// Create the W211 flux handler
	auto testFitFlux = make_shared<W211FitFluxHandler>();
	// Set the flux amplitude
	testFitFlux->setFluxAmplitude(1.0);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(network, surfacePos, grid);

	// Create a time
	double currTime = 1.0;

	// The array of concentration
	double newConcentration[5 * dof];

	// Initialize their values
	for (int i = 0; i < 5 * dof; i++) {
		newConcentration[i] = 0.0;
	}

	// The pointer to the grid point we want
	double *updatedConc = &newConcentration[0];
	double *updatedConcOffset = updatedConc + dof;

	// Update the concentrations at some grid points
	testFitFlux->computeIncidentFlux(currTime, updatedConcOffset, 1,
			surfacePos);
	updatedConcOffset = updatedConc + 2 * dof;
	testFitFlux->computeIncidentFlux(currTime, updatedConcOffset, 2,
			surfacePos);
	updatedConcOffset = updatedConc + 3 * dof;
	testFitFlux->computeIncidentFlux(currTime, updatedConcOffset, 3,
			surfacePos);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(newConcentration[9], 0.385429, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration[18], 0.286585, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration[27], 0.127985, 0.01);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()