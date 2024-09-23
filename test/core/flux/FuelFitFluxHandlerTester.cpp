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
	paramFile << "netParam=5 0 0 5 1" << std::endl
			  << "tempParam=1500" << std::endl;
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
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxXe = opts.getMaxImpurity();
	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		maxXe, maxV, maxI};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{maxXe + 1, maxV + 1, maxI + 1}};
	NetworkType network(maxSpeciesAmounts, subdivRatios, 1, opts);
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
	test::DOFView updatedConc("updatedConc", 5, dof);

	// Initialize their values
	Kokkos::parallel_for(
		Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {5, dof}),
		KOKKOS_LAMBDA(int i, int j) {
			conc(i, j) = 1.0e-5 * i;
			updatedConc(i, j) = 0.0;
		});

	// The pointer to the grid point we want
	auto concOffset = subview(conc, 1, Kokkos::ALL);
	auto updatedConcOffset = subview(updatedConc, 1, Kokkos::ALL);

	// Update the concentrations at some grid points
	testFitFlux->computeIncidentFlux(
		currTime, concOffset, updatedConcOffset, 1, surfacePos);
	concOffset = subview(conc, 2, Kokkos::ALL);
	updatedConcOffset = subview(updatedConc, 2, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(
		currTime, concOffset, updatedConcOffset, 2, surfacePos);
	concOffset = subview(conc, 3, Kokkos::ALL);
	updatedConcOffset = subview(updatedConc, 3, Kokkos::ALL);
	testFitFlux->computeIncidentFlux(
		currTime, concOffset, updatedConcOffset, 3, surfacePos);

	// Check the value at some grid points
	auto newConcentration =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConc);
	BOOST_REQUIRE_CLOSE(newConcentration(1, 1), 10000.0, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(2, 0), 10000.0, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration(3, 3), 0.25, 0.01);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
