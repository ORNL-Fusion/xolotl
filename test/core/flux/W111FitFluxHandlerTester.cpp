#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/flux/W111FitFluxHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/options/Options.h>
#include <xolotl/test/CommandLine.h>
#include <xolotl/util/MPIUtils.h>

using namespace std;
using namespace xolotl;
using namespace core;
using namespace flux;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * The test suite is responsible for testing the W111FitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE(W111FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkComputeIncidentFlux)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=9 0 0 0 0" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	util::mpiInit(cl.argc, cl.argv);
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 7; l++) {
		grid.push_back((double)l * 1.25);
	}
	// Specify the surface position
	int surfacePos = 0;

	// Create the network
	using NetworkType =
		network::PSIReactionNetwork<network::PSIFullSpeciesList>;
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe = opts.getMaxImpurity();
	NetworkType::AmountType maxD = opts.getMaxD();
	NetworkType::AmountType maxT = opts.getMaxT();
	NetworkType network({maxHe, maxD, maxT, maxV, maxI}, grid.size(), opts);
	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);
	// Get its size
	const int dof = network.getDOF();

	// Create the W111 flux handler
	auto testFitFlux = make_shared<W111FitFluxHandler>(opts);
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
	double* updatedConc = &newConcentration[0];
	double* updatedConcOffset = updatedConc + dof;

	// Update the concentrations at some grid points
	testFitFlux->computeIncidentFlux(
		currTime, updatedConcOffset, 1, surfacePos);
	updatedConcOffset = updatedConc + 2 * dof;
	testFitFlux->computeIncidentFlux(
		currTime, updatedConcOffset, 2, surfacePos);
	updatedConcOffset = updatedConc + 3 * dof;
	testFitFlux->computeIncidentFlux(
		currTime, updatedConcOffset, 3, surfacePos);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(newConcentration[9], 0.3168967, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration[18], 0.306857, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration[27], 0.1762458, 0.01);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
