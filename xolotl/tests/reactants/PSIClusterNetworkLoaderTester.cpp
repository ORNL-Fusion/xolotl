/*
 * PSIClusterNetworkLoaderTester.cpp
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSIClusterNetworkLoader.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <PSIClusterNetworkLoader.h>
#include <PSIClusterReactionNetwork.h>


using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the PSIClusterNetworkLoader. It
 * creates a string stream that contains each of the available PSICluster types
 * and checks that the loader returns a list with each type in it.
 */
BOOST_AUTO_TEST_SUITE(PSIClusterNetworkLoader_testSuite)

/** This operation checks the loader. */
BOOST_AUTO_TEST_CASE(checkLoading) {

	// Local Declarations
	shared_ptr<stringstream> networkStream(
			new stringstream(stringstream::in | stringstream::out));
	string singleHeString = "1 0 0 0.0 Infinity Infinity 8.269996 0.999 1.34\n";
	string singleVString =
			"0 50 0 Infinity 2.49000002 Infinity Infinity 0.888 2.345\n";
	string singleIString = "0 0 1 Infinity Infinity Infinity Infinity 0.7777 3.456\n";
	string mixedString =
			"1 50 0 6.160001 2.4900002 Infinity Infinity 6.789 4.5678\n";
	// This string is bad because it is one value short
	string badString = "1 2 3 4 5 6 7 8\n";
	bool caughtFlag = false;
	PSIClusterNetworkLoader loader = PSIClusterNetworkLoader();
	
	
	// Load the network stream. This simulates a file with single He, single
	// V, single I and one mixed-species cluster. They are mixed up here to test
	// the ability of the loader to order them.
	*networkStream << singleVString << mixedString << singleHeString
			<< singleIString;

	// Diagnostic information
	// @formatter:off
	BOOST_TEST_MESSAGE("CLUSTER DATA"
			<< "\nHe: "
			<< singleHeString
			<< "\nV: "
			<< singleVString
			<< "\nI: "
			<< singleIString
			<< "\n Mixed: "
			<< mixedString
			<< "\nFull Network data: \n"
			<< (*networkStream).str());
	// @formatter:off

	// Setup the Loader
	loader.setInputstream(networkStream);

	// Load the network
	shared_ptr<PSIClusterReactionNetwork> network = loader.load();

	// Check the network. It should not be empty
	auto props = network->getProperties();
	BOOST_TEST_MESSAGE("Number of properties = " << props.size());
	BOOST_REQUIRE(!props.empty());
	
	// Print the properties list to debug
	for (auto it = props.begin(); it != props.end(); it++) {
		printf("\"%s\" => \"%s\"\n", it->first.c_str(), it->second.c_str());
	}

	// Check the properties
	BOOST_TEST_MESSAGE("Maximum He Cluster Size = " << props["maxHeClusterSize"]);
	BOOST_REQUIRE(strtol(props["maxHeClusterSize"].c_str(),NULL,10) == 10);
	BOOST_TEST_MESSAGE("Maximum V Cluster Size = " << props["maxVClusterSize"]);
	BOOST_REQUIRE(strtol(props["maxVClusterSize"].c_str(),NULL,10) == 10);
	BOOST_TEST_MESSAGE("Maximum Interstitial Cluster Size = " << props["maxIClusterSize"]);
	BOOST_REQUIRE(strtol(props["maxIClusterSize"].c_str(),NULL,10) == 10);
	BOOST_TEST_MESSAGE("Maximum Mixed Species Cluster Size = " << props["maxMixedClusterSize"]);
	BOOST_REQUIRE(strtol(props["maxMixedClusterSize"].c_str(),NULL,10) == 10);
	BOOST_TEST_MESSAGE("Number of He clusters = " << props["numHeClusters"]);
	BOOST_REQUIRE(strtol(props["numHeClusters"].c_str(),NULL,10) == 1);
	BOOST_TEST_MESSAGE("Number of V clusters = " << props["numVClusters"]);
	BOOST_REQUIRE(strtol(props["numVClusters"].c_str(),NULL,10) == 1);
	BOOST_TEST_MESSAGE("Number of I clusters = " << props["numIClusters"]);
	BOOST_REQUIRE(strtol(props["numIClusters"].c_str(),NULL,10) == 1);

	// Check the reactants - He first
	std::shared_ptr<PSICluster> heCluster = static_pointer_cast<PSICluster>(network->get("He",1));
	BOOST_REQUIRE(heCluster->getSize() == 1);
	std::vector<double> bindingEnergies = heCluster->getBindingEnergies();
	BOOST_REQUIRE_CLOSE(bindingEnergies.at(0), 0.0, 0.00);
	BOOST_REQUIRE(bindingEnergies.at(1) == std::numeric_limits<double>::infinity());
	BOOST_REQUIRE(bindingEnergies.at(2) == std::numeric_limits<double>::infinity());
	BOOST_REQUIRE_CLOSE(bindingEnergies.at(3), 8.2699999, 0.0001);
	BOOST_REQUIRE_CLOSE(heCluster->getMigrationEnergy(),0.999,0.001);
	BOOST_REQUIRE_CLOSE(heCluster->getDiffusionFactor(),1.34,0.01);
	// V
	std::shared_ptr<PSICluster> vCluster = static_pointer_cast<PSICluster>(network->get("V",50));
	BOOST_REQUIRE(vCluster->getSize() == 50);
	bindingEnergies.clear();
	bindingEnergies = vCluster->getBindingEnergies();
	BOOST_REQUIRE(bindingEnergies[0] == std::numeric_limits<double>::infinity());
	BOOST_REQUIRE_CLOSE(bindingEnergies[1], 2.4900, 0.0001);
	BOOST_REQUIRE(bindingEnergies[2] == std::numeric_limits<double>::infinity());
	BOOST_REQUIRE(bindingEnergies[3] == std::numeric_limits<double>::infinity());
	BOOST_REQUIRE_CLOSE(vCluster->getMigrationEnergy(),0.888,0.001);
	BOOST_REQUIRE_CLOSE(vCluster->getDiffusionFactor(),2.345,0.001);
	// I
	std::shared_ptr<PSICluster> iCluster = static_pointer_cast<PSICluster>(network->get("I",1));
	BOOST_REQUIRE(iCluster->getSize() == 1);
	bindingEnergies.clear();
	bindingEnergies = iCluster->getBindingEnergies();
	BOOST_REQUIRE(bindingEnergies[0] == std::numeric_limits<double>::infinity());
	BOOST_REQUIRE(bindingEnergies[1] == std::numeric_limits<double>::infinity());
	BOOST_REQUIRE(bindingEnergies[2] == std::numeric_limits<double>::infinity());
	BOOST_REQUIRE(bindingEnergies[3] == std::numeric_limits<double>::infinity());
	BOOST_REQUIRE_CLOSE(iCluster->getMigrationEnergy(),0.7777,0.0001);
	BOOST_REQUIRE_CLOSE(iCluster->getDiffusionFactor(),3.456,0.001);
	// HeV
	std::vector<int> composition;
	composition.push_back(1);
	composition.push_back(50);
	composition.push_back(0);
	std::shared_ptr<PSICluster> heVCluster = static_pointer_cast<PSICluster>(network->getCompound("HeV",composition));
	BOOST_REQUIRE(heVCluster->getSize() == 51);
	bindingEnergies.clear();
	bindingEnergies = heVCluster->getBindingEnergies();
	BOOST_REQUIRE_CLOSE(bindingEnergies[0], 6.160001, 0.00001);
	BOOST_REQUIRE_CLOSE(bindingEnergies[1], 2.4900, 0.0001);
	BOOST_REQUIRE(bindingEnergies[2] == std::numeric_limits<double>::infinity());
	BOOST_REQUIRE(bindingEnergies[3] == std::numeric_limits<double>::infinity());
	BOOST_REQUIRE_CLOSE(heVCluster->getMigrationEnergy(),6.789,0.001);
	BOOST_REQUIRE_CLOSE(heVCluster->getDiffusionFactor(),4.5678,0.0001);

	// Reload the network stream with the bad string
	(*networkStream).clear();
	*networkStream << badString;
	// Make sure the exception is caught when loading the bad string
	try {
		loader.load();
	} catch (std::string error) {
		// Do nothing but flip the flag
		caughtFlag = true;
	}
	BOOST_REQUIRE(caughtFlag);

	return;
}
BOOST_AUTO_TEST_SUITE_END()
