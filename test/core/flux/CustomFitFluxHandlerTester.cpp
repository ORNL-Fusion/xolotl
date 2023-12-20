#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/flux/CustomFitFluxHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/options/Options.h>
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
 * The test suite is responsible for testing the CustomFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE(CustomFitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkComputeIncidentFlux)
{
	// Create the option to create a network
	xolotl::options::Options opts;

	// Create a file with flux profile data.
	std::ofstream fluxFile("tridyn.dat");
	fluxFile
		<< "He 1 0.07280658930900809" << std::endl
		<< "-0.007081476181737453 0.1261483878960484 0.046411385590726884 "
		   "-0.08026629198472447 0.039202798471250494 -0.011029329255965437 "
		   "0.002043027695798964 -0.00026258604985788386 "
		   "2.4021358408263736e-05 "
		   "-1.5804767473061154e-06 7.464764376419655e-08 "
		   "-2.49230642776284e-09 "
		   "5.6866990099240226e-11 -8.311086756674685e-13 "
		   "6.8300108211942295e-15 "
		   "-2.271813054793263e-17 25.053135779801487"
		<< std::endl
		<< "I 1 0.0019765160176424246" << std::endl
		<< "1.0443514769976874 6.961489658525669 -49.619639407212524 "
		   "307.35407166956384 -1385.4341206868723 3831.3312179118425 "
		   "-6796.230883904835 8173.956265416343 -6915.097674700663 "
		   "4197.699943058166 -1839.2348194439237 577.1583583299342 "
		   "-126.55005505770353 18.415889098836942 -1.598219349445743 "
		   "0.06260557665223535 3.3074402233784164"
		<< std::endl;
	fluxFile.close();

	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 0 0 10 6" << std::endl
			  << "fluxDepthProfileFilePath=tridyn.dat" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	util::mpiInit(cl.argc, cl.argv);
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 6; l++) {
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
	// Get its size
	const int dof = network.getDOF();

	// Create the Custom flux handler
	auto testFitFlux = make_shared<CustomFitFluxHandler>(opts);
	// Set the flux amplitude
	testFitFlux->setFluxAmplitude(1.0);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(network, surfacePos, grid);

	// Check the cluster Ids
	auto idVector = testFitFlux->getFluxIndices();
	BOOST_REQUIRE_EQUAL(idVector.size(), 2);
	BOOST_REQUIRE_EQUAL(idVector[0], 16);
	BOOST_REQUIRE_EQUAL(idVector[1], 0);

	// Create a time
	double currTime = 1.0;

	// The array of concentration
	test::DOFView conc("conc", 6, dof);

	// The pointer to the grid point we want
	auto updatedConcOffset = subview(conc, 1, Kokkos::ALL);

	// Update the concentrations at some grid points
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 1, surfacePos);
	updatedConcOffset = subview(conc, 2, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 2, surfacePos);

	// Check the value at some grid points
	auto newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, conc);
	BOOST_REQUIRE_CLOSE(newConcentration(1, 0), 0.0015776, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(1, 16), 0.01840469, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(1, 46), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(2, 0), 3.60006e-06, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(2, 16), 0.0398406, 0.01);

	return;
}

BOOST_AUTO_TEST_CASE(checkFluence)
{
	// Create the option to create a network
	xolotl::options::Options opts;

	// Create a file with flux profile data.
	std::ofstream fluxFile("tridyn.dat");
	fluxFile
		<< "He 1 0.07280658930900809" << std::endl
		<< "-0.007081476181737453 0.1261483878960484 0.046411385590726884 "
		   "-0.08026629198472447 0.039202798471250494 -0.011029329255965437 "
		   "0.002043027695798964 -0.00026258604985788386 "
		   "2.4021358408263736e-05 "
		   "-1.5804767473061154e-06 7.464764376419655e-08 "
		   "-2.49230642776284e-09 "
		   "5.6866990099240226e-11 -8.311086756674685e-13 "
		   "6.8300108211942295e-15 "
		   "-2.271813054793263e-17 25.053135779801487"
		<< std::endl
		<< "I 1 0.0019765160176424246" << std::endl
		<< "1.0443514769976874 6.961489658525669 -49.619639407212524 "
		   "307.35407166956384 -1385.4341206868723 3831.3312179118425 "
		   "-6796.230883904835 8173.956265416343 -6915.097674700663 "
		   "4197.699943058166 -1839.2348194439237 577.1583583299342 "
		   "-126.55005505770353 18.415889098836942 -1.598219349445743 "
		   "0.06260557665223535 3.3074402233784164"
		<< std::endl;
	fluxFile.close();

	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 0 0 10 6" << std::endl
			  << "fluxDepthProfileFilePath=tridyn.dat" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
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
	// Get its size
	const int dof = network.getDOF();

	// Create the W100 flux handler
	auto testFitFlux = make_shared<CustomFitFluxHandler>(opts);
	// Set the flux amplitude
	testFitFlux->setFluxAmplitude(1.0);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(network, surfacePos, grid);

	// Check that the fluence is 0.0 at the beginning
	auto fluence = testFitFlux->getFluence();
	BOOST_REQUIRE_EQUAL(fluence[0], 0.0);

	// Increment the fluence
	testFitFlux->incrementFluence(1.0e-8);

	// Check that the fluence is not 0.0 anymore
	fluence = testFitFlux->getFluence();
	BOOST_REQUIRE_EQUAL(fluence[0], 1.0e-8);
	BOOST_REQUIRE_CLOSE(fluence[1], 7.28065893e-10, 0.001);
	BOOST_REQUIRE_CLOSE(fluence[2], 1.976516e-11, 0.001);

	// Check setFluence
	std::vector<double> f = {1.0, 5.0, 0.2};
	testFitFlux->setFluence(f);
	fluence = testFitFlux->getFluence();
	BOOST_REQUIRE_EQUAL(fluence[0], f[0]);
	BOOST_REQUIRE_EQUAL(fluence[1], f[1]);
	BOOST_REQUIRE_EQUAL(fluence[2], f[2]);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
