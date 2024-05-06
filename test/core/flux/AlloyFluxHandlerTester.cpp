#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/flux/A800H5MeVFluxHandler.h>
#include <xolotl/core/flux/A800H9MeVFluxHandler.h>
#include <xolotl/core/flux/A800HNeutronFluxHandler.h>
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
 * The test suite is responsible for testing the A800H*FluxHandler.
 */
BOOST_AUTO_TEST_SUITE(A800HFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(check800H5MeV)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=50 50 0 6 4" << std::endl;
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
	using NetworkType = network::AlloyReactionNetwork;
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxSize = opts.getMaxImpurity();
	NetworkType network(
		{maxV, maxSize, maxSize, maxI, maxSize, maxSize}, 0, opts);
	// Get its size
	const int dof = network.getDOF();

	// Create the iron flux handler
	auto testFlux = make_shared<A800H5MeVFluxHandler>(opts);
	// Set the flux amplitude
	testFlux->setFluxAmplitude(1.0);
	// Initialize the flux handler
	testFlux->initializeFluxHandler(network, surfacePos, grid);

	// Create a time
	double currTime = 1.0;

	// The array of concentration
	test::DOFView conc("conc", 7, dof);

	// The pointer to the grid point we want
	auto updatedConcOffset = subview(conc, 1, Kokkos::ALL);

	// Update the concentrations
	testFlux->computeIncidentFlux(currTime, updatedConcOffset, 1, surfacePos);

	// Check the value at some grid points
	auto newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(newConcentration[93], 4.49925, 0.01); // I_1
	BOOST_REQUIRE_CLOSE(newConcentration[96], 0.2121075, 0.01); // I_4
	BOOST_REQUIRE_CLOSE(newConcentration[47], 0.071988, 0.01); // Perfect_5
	BOOST_REQUIRE_CLOSE(newConcentration[87], 0.0, 0.01); // Perfect_45
	BOOST_REQUIRE_CLOSE(newConcentration[0], 0.017997, 0.01); // Frank_5
	BOOST_REQUIRE_CLOSE(newConcentration[40], 0.0, 0.01); // Frank_45
	BOOST_REQUIRE_CLOSE(newConcentration[185], 4.49925, 0.01); // V_1
	BOOST_REQUIRE_CLOSE(newConcentration[190], 0.0835575, 0.01); // V_6
	BOOST_REQUIRE_CLOSE(newConcentration[97], 0.034892, 0.01); // Faulted_7
	BOOST_REQUIRE_CLOSE(newConcentration[135], 0.0, 0.01); // Faulted_45

	updatedConcOffset = subview(conc, 4, Kokkos::ALL);

	// Update the concentrations
	testFlux->computeIncidentFlux(currTime, updatedConcOffset, 4, surfacePos);

	// Check the value at some grid points
	deep_copy(newConcentration, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(newConcentration[93], 0.0, 0.01); // I_1
	BOOST_REQUIRE_CLOSE(newConcentration[96], 0.0, 0.01); // I_4
	BOOST_REQUIRE_CLOSE(newConcentration[47], 0.0, 0.01); // Perfect_5
	BOOST_REQUIRE_CLOSE(newConcentration[87], 0.0, 0.01); // Perfect_45
	BOOST_REQUIRE_CLOSE(newConcentration[0], 0.0, 0.01); // Frank_5
	BOOST_REQUIRE_CLOSE(newConcentration[40], 0.0, 0.01); // Frank_45
	BOOST_REQUIRE_CLOSE(newConcentration[185], 0.0, 0.01); // V_1
	BOOST_REQUIRE_CLOSE(newConcentration[190], 0.0, 0.01); // V_6
	BOOST_REQUIRE_CLOSE(newConcentration[97], 0.0, 0.01); // Faulted_7
	BOOST_REQUIRE_CLOSE(newConcentration[135], 0.0, 0.01); // Faulted_45

	return;
}

BOOST_AUTO_TEST_CASE(check800H9MeV)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=50 50 0 6 4" << std::endl;
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
	using NetworkType = network::AlloyReactionNetwork;
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxSize = opts.getMaxImpurity();
	NetworkType network(
		{maxV, maxSize, maxSize, maxI, maxSize, maxSize}, 0, opts);
	// Get its size
	const int dof = network.getDOF();

	// Create the iron flux handler
	auto testFlux = make_shared<A800H9MeVFluxHandler>(opts);
	// Set the flux amplitude
	testFlux->setFluxAmplitude(1.0);
	// Initialize the flux handler
	testFlux->initializeFluxHandler(network, surfacePos, grid);

	// Create a time
	double currTime = 1.0;

	// The array of concentration
	test::DOFView conc("conc", 7, dof);

	// The pointer to the grid point we want
	auto updatedConcOffset = subview(conc, 1, Kokkos::ALL);

	// Update the concentrations
	testFlux->computeIncidentFlux(currTime, updatedConcOffset, 1, surfacePos);

	// Check the value at some grid points
	auto newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(newConcentration[93], 4.49925, 0.01); // I_1
	BOOST_REQUIRE_CLOSE(newConcentration[96], 0.2121075, 0.01); // I_4
	BOOST_REQUIRE_CLOSE(newConcentration[47], 0.071988, 0.01); // Perfect_5
	BOOST_REQUIRE_CLOSE(newConcentration[87], 0.0, 0.01); // Perfect_45
	BOOST_REQUIRE_CLOSE(newConcentration[0], 0.017997, 0.01); // Frank_5
	BOOST_REQUIRE_CLOSE(newConcentration[40], 0.0, 0.01); // Frank_45
	BOOST_REQUIRE_CLOSE(newConcentration[185], 4.49925, 0.01); // V_1
	BOOST_REQUIRE_CLOSE(newConcentration[190], 0.0835575, 0.01); // V_6
	BOOST_REQUIRE_CLOSE(newConcentration[97], 0.034892, 0.01); // Faulted_7
	BOOST_REQUIRE_CLOSE(newConcentration[135], 0.0, 0.01); // Faulted_45

	updatedConcOffset = subview(conc, 4, Kokkos::ALL);

	// Update the concentrations
	testFlux->computeIncidentFlux(currTime, updatedConcOffset, 4, surfacePos);

	// Check the value at some grid points
	deep_copy(newConcentration, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(newConcentration[93], 0.0, 0.01); // I_1
	BOOST_REQUIRE_CLOSE(newConcentration[96], 0.0, 0.01); // I_4
	BOOST_REQUIRE_CLOSE(newConcentration[47], 0.0, 0.01); // Perfect_5
	BOOST_REQUIRE_CLOSE(newConcentration[87], 0.0, 0.01); // Perfect_45
	BOOST_REQUIRE_CLOSE(newConcentration[0], 0.0, 0.01); // Frank_5
	BOOST_REQUIRE_CLOSE(newConcentration[40], 0.0, 0.01); // Frank_45
	BOOST_REQUIRE_CLOSE(newConcentration[185], 0.0, 0.01); // V_1
	BOOST_REQUIRE_CLOSE(newConcentration[190], 0.0, 0.01); // V_6
	BOOST_REQUIRE_CLOSE(newConcentration[97], 0.0, 0.01); // Faulted_7
	BOOST_REQUIRE_CLOSE(newConcentration[135], 0.0, 0.01); // Faulted_45

	return;
}

BOOST_AUTO_TEST_CASE(check800HNeutron)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=50 50 0 6 4" << std::endl;
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
	using NetworkType = network::AlloyReactionNetwork;
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxSize = opts.getMaxImpurity();
	NetworkType network(
		{maxV, maxSize, maxSize, maxI, maxSize, maxSize}, 0, opts);
	// Get its size
	const int dof = network.getDOF();

	// Create the iron flux handler
	auto testFlux = make_shared<A800H9MeVFluxHandler>(opts);
	// Set the flux amplitude
	testFlux->setFluxAmplitude(1.0);
	// Initialize the flux handler
	testFlux->initializeFluxHandler(network, surfacePos, grid);

	// Create a time
	double currTime = 1.0;

	// The array of concentration
	test::DOFView conc("conc", 7, dof);

	// The pointer to the grid point we want
	auto updatedConcOffset = subview(conc, 1, Kokkos::ALL);

	// Update the concentrations
	testFlux->computeIncidentFlux(currTime, updatedConcOffset, 1, surfacePos);

	// Check the value at some grid points
	auto newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(newConcentration[93], 4.49925, 0.01); // I_1
	BOOST_REQUIRE_CLOSE(newConcentration[96], 0.2121075, 0.01); // I_4
	BOOST_REQUIRE_CLOSE(newConcentration[47], 0.071988, 0.01); // Perfect_5
	BOOST_REQUIRE_CLOSE(newConcentration[87], 0.0, 0.01); // Perfect_45
	BOOST_REQUIRE_CLOSE(newConcentration[0], 0.017997, 0.01); // Frank_5
	BOOST_REQUIRE_CLOSE(newConcentration[40], 0.0, 0.01); // Frank_45
	BOOST_REQUIRE_CLOSE(newConcentration[185], 4.49925, 0.01); // V_1
	BOOST_REQUIRE_CLOSE(newConcentration[190], 0.0835575, 0.01); // V_6
	BOOST_REQUIRE_CLOSE(newConcentration[97], 0.034892, 0.01); // Faulted_7
	BOOST_REQUIRE_CLOSE(newConcentration[135], 0.0, 0.01); // Faulted_45

	updatedConcOffset = subview(conc, 4, Kokkos::ALL);

	// Update the concentrations
	testFlux->computeIncidentFlux(currTime, updatedConcOffset, 4, surfacePos);

	// Check the value at some grid points
	deep_copy(newConcentration, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(newConcentration[93], 0.0, 0.01); // I_1
	BOOST_REQUIRE_CLOSE(newConcentration[96], 0.0, 0.01); // I_4
	BOOST_REQUIRE_CLOSE(newConcentration[47], 0.0, 0.01); // Perfect_5
	BOOST_REQUIRE_CLOSE(newConcentration[87], 0.0, 0.01); // Perfect_45
	BOOST_REQUIRE_CLOSE(newConcentration[0], 0.0, 0.01); // Frank_5
	BOOST_REQUIRE_CLOSE(newConcentration[40], 0.0, 0.01); // Frank_45
	BOOST_REQUIRE_CLOSE(newConcentration[185], 0.0, 0.01); // V_1
	BOOST_REQUIRE_CLOSE(newConcentration[190], 0.0, 0.01); // V_6
	BOOST_REQUIRE_CLOSE(newConcentration[97], 0.0, 0.01); // Faulted_7
	BOOST_REQUIRE_CLOSE(newConcentration[135], 0.0, 0.01); // Faulted_45

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
