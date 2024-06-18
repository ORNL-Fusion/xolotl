#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/flux/FuelFitFluxHandler.h>
#include <xolotl/options/ConfOptions.h>
#include <xolotl/test/CommandLine.h>
#include <xolotl/test/Util.h>
#include <xolotl/util/MPIUtils.h>

using namespace std;
using namespace xolotl;
using namespace core;
using namespace flux;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE(FuelFitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkComputeIncidentFlux)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=100 0 0 0 0" << std::endl;
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
	using NetworkType = network::NEReactionNetwork;
	NetworkType::AmountType maxXe = opts.getMaxImpurity();
	NetworkType network({maxXe}, grid.size(), opts);
	// Get its size
	const int dof = network.getDOF();

	// Create the fuel flux handler
	auto testFitFlux = make_shared<FuelFitFluxHandler>(opts);
	// Set the flux amplitude
	testFitFlux->setFluxAmplitude(1.0);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(network, surfacePos, grid);

	// Create a time
	double currTime = 1.0;

	// The array of concentration
	test::DOFView conc("conc", 5, dof);

	// The pointer to the grid point we want
	auto updatedConcOffset = subview(conc, 1, Kokkos::ALL);

	// Update the concentrations at some grid points
	testFitFlux->computeIncidentFlux(
		currTime, updatedConcOffset, 1, surfacePos);
	updatedConcOffset = subview(conc, 2, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(
		currTime, updatedConcOffset, 2, surfacePos);
	updatedConcOffset = subview(conc, 3, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(
		currTime, updatedConcOffset, 3, surfacePos);

	// Check the value at some grid points
	auto newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, conc);
	BOOST_REQUIRE_CLOSE(newConcentration(1, 0), 1.0, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(2, 0), 1.0, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(3, 0), 1.0, 0.01);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
