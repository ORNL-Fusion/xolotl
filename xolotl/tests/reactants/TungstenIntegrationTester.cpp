#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSIClusterNetworkLoader.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <PSIClusterNetworkLoader.h>
#include <PSIClusterReactionNetwork.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>

using namespace std;
using namespace xolotlCore;

/**
 * The test suite configuration
 */BOOST_AUTO_TEST_SUITE (TungstenIntegrationTester_testSuite)

/**
 * I just noticed that there is no unit test for checking the partial
 * derivatives computed by Xolotl. I think I did this because there is no good
 * way to unit test them other than computing them and then sticking that array
 * in for a direct comparison. Ideally there would be some analytic version to
 * which we could compare.
 *
 * The simplest way to test these is to use PETSc’s -snes_type test option and
 * grep/diff the output.
 *
 * I just modified the test below to checks the convenience method I added to
 * get the partials in a pre-allocated array against the original. They should
 * give the same results.
 *
 * ~JJB 20140525 10:19
 *
 */

/**
 * This operation checks the fluxs from the reactant as best as is possible
 * given that it requires external data.
 */
BOOST_AUTO_TEST_CASE(checkGetReactantFluxesAndParials) {
	// Local Declarations
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/tungsten.txt");
	string networkFilename = sourceDir + pathToFile;

	BOOST_TEST_MESSAGE(
			"TungstenIntegrationTester Message: Network filename is: " << networkFilename);

	// Load the input file from the master task
	shared_ptr<istream> networkStream = make_shared < ifstream
			> (networkFilename);

	// Create a network loader and set the istream on every MPI task
	shared_ptr<PSIClusterNetworkLoader> networkLoader = make_shared
			< PSIClusterNetworkLoader
			> (std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	networkLoader->setInputstream(networkStream);
	// Load the network
	shared_ptr<ReactionNetwork> network = networkLoader->load();

	BOOST_TEST_MESSAGE("TungstenIntegrationTester Message: Network loaded");

	// Get everything from the network and set the temperature
	double temperature = 1000.0;
	int nReactants = network->size();
	auto reactants = network->getAll();
	network->setTemperature(temperature);

	BOOST_TEST_MESSAGE("TungstenIntegrationTester Message: "
			<< "Size of the network is: "
			<< nReactants);

	// Create an array to test the second partial derivative routine.
	auto secondPartials = std::vector<double>(nReactants,0.0);

	BOOST_TEST_MESSAGE("Check partial derivatives.");
	for (int i = 0; i < nReactants; ++i) {
		auto reactant = (PSICluster *) reactants->at(i);
		double flux = reactant->getTotalFlux(temperature);
		// Get the partials using method 1
		auto partials = reactant->getPartialDerivatives(temperature);
		// Get the partials using method 2
		reactant->getPartialDerivatives(temperature,secondPartials);
		// Compare the two arrays of partial derivatives
		for (int j = 0; j < nReactants; ++j) {
			BOOST_REQUIRE_CLOSE(partials[j],secondPartials[j],1.0);
		}
		// Zero the partials array
		std::fill(secondPartials.begin(),secondPartials.end(),0.0);
	}

}
BOOST_AUTO_TEST_SUITE_END()
