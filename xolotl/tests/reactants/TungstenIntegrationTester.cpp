#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSIClusterNetworkLoader.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <PSIClusterNetworkLoader.h>
#include <PSIClusterReactionNetwork.h>
#include "../../XolotlConfig.h"

using namespace std;
using namespace xolotlCore;

/**
 * The test suite configuration
 */
BOOST_AUTO_TEST_SUITE (TungstenIntegrationTester_testSuite)

/**
 * This operation checks the fluxs from the reactant as best as is possible
 * given that it requires external data.
 */
BOOST_AUTO_TEST_CASE(checkGetReactantFluxes) {
	// Local Declarations
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/tungsten.txt");
	string networkFilename = sourceDir + pathToFile;

	std::cout << networkFilename << std::endl;
	// Load the input file from the master task
	shared_ptr<std::istream> networkStream;
	networkStream.reset(new std::ifstream(networkFilename));

	// Create a network loader and set the istream on every MPI task
	std::shared_ptr<PSIClusterNetworkLoader> networkLoader(
			new PSIClusterNetworkLoader());
	networkLoader->setInputstream(networkStream);
	// Load the network
	shared_ptr<ReactionNetwork> network = networkLoader->load();

	cout << "Network loaded." << endl;

	int nReactants = network->size();
	auto reactants = network->getAll();
	cout << "Size of the network is " << nReactants << "\n";
	for (int i = 0; i < 9; ++i) {
		shared_ptr<PSICluster> reactant = dynamic_pointer_cast < PSICluster > (reactants->at(i));
		double flux = reactant->getTotalFlux(273.0);
	}

}

BOOST_AUTO_TEST_SUITE_END()
