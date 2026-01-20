#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/flux/PulsedFitFluxHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
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
 * The test suite is responsible for testing the PulsedFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE(PulsedFitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkComputeIncidentFlux)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=0 0 0 2 2" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	util::mpiInit(cl.argc, cl.argv);
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 100; l++) {
		grid.push_back((double)l * 50.0);
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
	// Get its size
	const int dof = network.getDOF();

	// Create the ulsed flux handler
	auto testFitFlux = make_shared<PulsedFitFluxHandler>(opts);
	// Set the flux amplitude and pulse parameters
	testFitFlux->setFluxAmplitude(1.0);
	testFitFlux->setPulseTime(1.0e-3);
	testFitFlux->setProportion(0.2);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(network, surfacePos, grid);

	// Create a time
	double currTime = 1.0e-4;

	// The array of concentration
	test::DOFView conc("conc", 100, dof);

	// The pointer to the grid point we want
	auto updatedConcOffset = subview(conc, 50, Kokkos::ALL);

	// Update the concentrations at some grid point
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 50, surfacePos);

	// Check the value at some grid point
	auto newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(newConcentration[0], 5.0295072885924443e-08, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration[2], 5.0295072885924443e-08, 0.01);

	// Check that the flux is 0.0 at later time
	currTime = 8.0e-4;
	updatedConcOffset = subview(conc, 22, Kokkos::ALL);

	// Update the concentrations at some grid point
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 22, surfacePos);
	newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);

	// Check the value at some grid point
	BOOST_REQUIRE_CLOSE(newConcentration[0], 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration[2], 0.0, 0.01);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
