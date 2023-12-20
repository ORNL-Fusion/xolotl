#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/flux/W100FitFluxHandler.h>
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
 * The test suite is responsible for testing the W100FitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE(W100FitFluxHandlerTester_testSuite)

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
	// Get its size
	const int dof = network.getDOF();

	// Create the W100 flux handler
	auto testFitFlux = make_shared<W100FitFluxHandler>(opts);
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
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 1, surfacePos);
	updatedConcOffset = subview(conc, 2, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 2, surfacePos);
	updatedConcOffset = subview(conc, 3, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 3, surfacePos);

	// Check the value at some grid points
	auto newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, conc);
	BOOST_REQUIRE_CLOSE(newConcentration(1, 0), 0.444777, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(2, 0), 0.247638, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(3, 0), 0.10758, 0.01);

	return;
}

BOOST_AUTO_TEST_CASE(checkComputeIncidentFluxNoGrid)
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
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	// Create a grid
	std::vector<double> grid;
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
	auto testFitFlux = make_shared<W100FitFluxHandler>(opts);
	// Set the flux amplitude
	testFitFlux->setFluxAmplitude(1.0);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(network, surfacePos, grid);

	// Create a time
	double currTime = 1.0;

	// The array of concentration
	Kokkos::View<double**> conc("conc", 1, dof);

	// The pointer to the grid point we want
	auto updatedConcOffset = subview(conc, 0, Kokkos::ALL);

	// Update the concentrations at some grid points
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 0, surfacePos);

	// Check the value at some grid points
	auto newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(newConcentration[0], 1.0, 0.01);

	return;
}

BOOST_AUTO_TEST_CASE(checkFluence)
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
	auto testFitFlux = make_shared<W100FitFluxHandler>(opts);
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

	return;
}

BOOST_AUTO_TEST_CASE(checkFluxAmplitude)
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
	auto testFitFlux = make_shared<W100FitFluxHandler>(opts);

	// Set the factor to change the flux amplitude
	testFitFlux->setFluxAmplitude(2.5);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(network, surfacePos, grid);

	// Check the value of the flux amplitude
	BOOST_REQUIRE_EQUAL(testFitFlux->getFluxAmplitude(), 2.5);

	// Create a time
	double currTime = 1.0;

	// The array of concentration
	test::DOFView conc("conc", 5, dof);

	// The pointer to the grid point we want
	auto updatedConcOffset = subview(conc, 1, Kokkos::ALL);

	// Update the concentrations at some grid points
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 1, surfacePos);
	updatedConcOffset = subview(conc, 2, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 2, surfacePos);
	updatedConcOffset = subview(conc, 3, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 3, surfacePos);

	// Check the value at some grid points
	auto newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, conc);
	BOOST_REQUIRE_CLOSE(newConcentration(1, 0), 1.111943, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(2, 0), 0.619095, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(3, 0), 0.268961, 0.01);

	return;
}

BOOST_AUTO_TEST_CASE(checkTimeProfileFlux)
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

	// Create a file with a time profile for the flux
	// First column with the time and the second with
	// the amplitude (in He/nm2/s) at that time.
	std::string fluxFile = "fluxFile.dat";
	std::ofstream writeFluxFile(fluxFile);
	writeFluxFile << "0.0 1000.0 \n"
					 "1.0 4000.0 \n"
					 "2.0 2000.0 \n"
					 "3.0 3000.0 \n"
					 "4.0 0.0";
	writeFluxFile.close();

	auto testFitFlux = make_shared<W100FitFluxHandler>(opts);
	// Initialize the time profile for the flux handler
	testFitFlux->initializeTimeProfile(fluxFile);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(network, surfacePos, grid);

	// Check flux indices
	auto indices = testFitFlux->getFluxIndices();
	BOOST_REQUIRE_EQUAL(indices.size(), 1);
	BOOST_REQUIRE_EQUAL(indices[0], 0);

	// Create a time
	double currTime = 0.5;

	// Check the instant flux
	auto instantFlux = testFitFlux->getInstantFlux(currTime);
	BOOST_REQUIRE_EQUAL(instantFlux.size(), 1);
	BOOST_REQUIRE_CLOSE(instantFlux[0], 2500.0, 0.01);

	// The array of concentration
	test::DOFView conc("conc", 5, dof);

	// The pointer to the grid point we want
	auto updatedConcOffset = subview(conc, 1, Kokkos::ALL);

	// Update the concentrations at some grid points
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 1, surfacePos);
	updatedConcOffset = subview(conc, 2, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 2, surfacePos);
	updatedConcOffset = subview(conc, 3, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 3, surfacePos);

	// Check the value at some grid points
	auto newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, conc);
	BOOST_REQUIRE_CLOSE(newConcentration(1, 0), 1111.94, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(2, 0), 619.095, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(3, 0), 268.961, 0.01);
	// Check the value of the flux amplitude
	BOOST_REQUIRE_EQUAL(testFitFlux->getFluxAmplitude(), 2500.0);

	// Change the current time
	currTime = 3.5;

	// Check the instant flux
	instantFlux = testFitFlux->getInstantFlux(currTime);
	BOOST_REQUIRE_EQUAL(instantFlux.size(), 1);
	BOOST_REQUIRE_CLOSE(instantFlux[0], 1500.0, 0.01);

	// Reinitialize their values
	conc = test::DOFView("conc", 5, dof);

	// Update the concentrations at some grid points
	updatedConcOffset = subview(conc, 1, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 1, surfacePos);
	updatedConcOffset = subview(conc, 2, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 2, surfacePos);
	updatedConcOffset = subview(conc, 3, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(currTime, Kokkos::View<const double*>(),
		updatedConcOffset, 3, surfacePos);

	// Check the value at some grid points
	deep_copy(newConcentration, conc);
	BOOST_REQUIRE_CLOSE(newConcentration(1, 0), 667.166, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(2, 0), 371.457, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(3, 0), 161.377, 0.01);
	// Check the value of the flux amplitude
	BOOST_REQUIRE_EQUAL(testFitFlux->getFluxAmplitude(), 1500.0);

	// Remove the created file
	std::remove(fluxFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
