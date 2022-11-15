#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/flux/AlloyFitFluxHandler.h>
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
 * The test suite is responsible for testing the AlloyFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE(AlloyFitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkComputeIncidentFlux)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=50 0 0 6 4" << std::endl;
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
	auto testFitFlux = make_shared<AlloyFitFluxHandler>(opts);
	// Set the flux amplitude
	testFitFlux->setFluxAmplitude(1.0);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(network, surfacePos, grid);

	// Create a time
	double currTime = 1.0;

	// The array of concentration
	double newConcentration[7 * dof];

	// Initialize their values
	for (int i = 0; i < 7 * dof; i++) {
		newConcentration[i] = 0.0;
	}

	// The pointer to the grid point we want
	double* updatedConc = &newConcentration[0];
	double* updatedConcOffset = updatedConc + dof;

	// Update the concentrations
	testFitFlux->computeIncidentFlux(
		currTime, nullptr, updatedConcOffset, 1, surfacePos);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(newConcentration[88 + dof], 4.49925, 0.01); // I_1
	BOOST_REQUIRE_CLOSE(newConcentration[91 + dof], 0.2121075, 0.01); // I_4
	BOOST_REQUIRE_CLOSE(
		newConcentration[47 + dof], 0.071988, 0.01); // Perfect_5
	BOOST_REQUIRE_CLOSE(newConcentration[87 + dof], 0.0, 0.01); // Perfect_45
	BOOST_REQUIRE_CLOSE(newConcentration[dof], 0.017997, 0.01); // Frank_5
	BOOST_REQUIRE_CLOSE(newConcentration[40 + dof], 0.0, 0.01); // Frank_45
	BOOST_REQUIRE_CLOSE(newConcentration[180 + dof], 4.49925, 0.01); // V_1
	BOOST_REQUIRE_CLOSE(newConcentration[185 + dof], 0.0835575, 0.01); // V_6
	BOOST_REQUIRE_CLOSE(
		newConcentration[92 + dof], 0.034892, 0.01); // Faulted_7
	BOOST_REQUIRE_CLOSE(newConcentration[130 + dof], 0.0, 0.01); // Faulted_45

	updatedConcOffset = updatedConc + 4 * dof;

	// Update the concentrations
	testFitFlux->computeIncidentFlux(
		currTime, nullptr, updatedConcOffset, 4, surfacePos);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(newConcentration[88 + 4 * dof], 0.0, 0.01); // I_1
	BOOST_REQUIRE_CLOSE(newConcentration[91 + 4 * dof], 0.0, 0.01); // I_4
	BOOST_REQUIRE_CLOSE(newConcentration[47 + 4 * dof], 0.0, 0.01); // Perfect_5
	BOOST_REQUIRE_CLOSE(
		newConcentration[87 + 4 * dof], 0.0, 0.01); // Perfect_45
	BOOST_REQUIRE_CLOSE(newConcentration[4 * dof], 0.0, 0.01); // Frank_5
	BOOST_REQUIRE_CLOSE(newConcentration[40 + 4 * dof], 0.0, 0.01); // Frank_45
	BOOST_REQUIRE_CLOSE(newConcentration[180 + 4 * dof], 0.0, 0.01); // V_1
	BOOST_REQUIRE_CLOSE(newConcentration[185 + 4 * dof], 0.0, 0.01); // V_6
	BOOST_REQUIRE_CLOSE(newConcentration[92 + 4 * dof], 0.0, 0.01); // Faulted_7
	BOOST_REQUIRE_CLOSE(
		newConcentration[130 + 4 * dof], 0.0, 0.01); // Faulted_45

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
