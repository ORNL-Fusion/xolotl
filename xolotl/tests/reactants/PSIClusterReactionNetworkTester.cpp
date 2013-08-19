/*
 * PSIClusterTester.cpp
 *
 *  Created on: May 6, 2013
 *      Author: Jay Jay Billings
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include <map>
#include <memory>
#include <typeinfo>
#include <limits>
#include <math.h>
#include "SimpleReactionNetwork.h"
#include <HeVCluster.h>
#include <HeCluster.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the ReactionNetwork
 */
BOOST_AUTO_TEST_SUITE(ReactionNetwork_testSuite)

/**
 * This operation tests the copy constructor.
 */
BOOST_AUTO_TEST_CASE(checkCopying) {
	
	PSIClusterReactionNetwork network;
	
	// Set some properties
	(*network.properties)["numMixedClusters"] = "4";
	
	// Add a reactant
	network.reactants->push_back(std::shared_ptr<Reactant>(new Reactant));
	network.reactants->at(0)->setConcentration(50.0);
	
	
	// Copy the network
	PSIClusterReactionNetwork network2 = network;
	
	// Check that the ReactionNetwork fields are copied
	BOOST_REQUIRE_NE(network.properties.get(), network2.properties.get());
	BOOST_REQUIRE_EQUAL(network.properties->size(), network2.properties->size());
	
	BOOST_REQUIRE_NE(network.reactants.get(), network2.reactants.get());
	BOOST_REQUIRE_EQUAL(network.reactants->size(), network2.reactants->size());
	
	// Change the values of the copied ReactionNetwork
	// This shouldn't modify the Reactants contained inside the
	// first ReactionNetwork.
	
	(*network2.properties)["numMixedClusters"] = "5";
	network2.reactants->at(0)->setConcentration(52.0);
	
	BOOST_REQUIRE_NE(network.properties->at("numMixedClusters"),
		network2.properties->at("numMixedClusters"));
	BOOST_REQUIRE_CLOSE(network.reactants->at(0)->getConcentration(), 50.0, 1e-5);
	BOOST_REQUIRE_CLOSE(network2.reactants->at(0)->getConcentration(), 52.0, 1e-5);
}


BOOST_AUTO_TEST_CASE(clusterIndexConversions) {
	
	shared_ptr<ReactionNetwork> network = testUtils::getSimpleReactionNetwork();
	
	// Convert each index to a cluster map and back
	
	printf("i\tHe\tV\tI\n");
	printf("===\t===\t===\t===\n");
	
	int reactionLength = network->reactants->size();
	
	for (int index = 0; index < reactionLength; index++) {
		std::map<std::string, int> speciesMap = network->toClusterMap(index);
		printf("%d\t%d\t%d\t%d\n", index,
			speciesMap["He"], speciesMap["V"], speciesMap["I"]);
		int convertedIndex = network->toClusterIndex(speciesMap);
		
		BOOST_REQUIRE_EQUAL(convertedIndex, index);
	}
}

BOOST_AUTO_TEST_CASE(toClusterMap) {
	shared_ptr<ReactionNetwork> network = testUtils::getSimpleReactionNetwork();
	
	std::map<std::string, int> cluster;
	
	// Test a couple of the HeClusters
	
	cluster = network->toClusterMap(0);
	BOOST_REQUIRE_EQUAL(cluster["He"], 1);
	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
	
	cluster = network->toClusterMap(9);
	BOOST_REQUIRE_EQUAL(cluster["He"], 10);
	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
	
	// Test VClusters
	
	cluster = network->toClusterMap(10);
	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
	BOOST_REQUIRE_EQUAL(cluster["V"], 1);
	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
	
	cluster = network->toClusterMap(19);
	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
	BOOST_REQUIRE_EQUAL(cluster["V"], 10);
	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
	
	// Test IClusters
	
	cluster = network->toClusterMap(20);
	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
	BOOST_REQUIRE_EQUAL(cluster["I"], 1);
	
	cluster = network->toClusterMap(29);
	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
	BOOST_REQUIRE_EQUAL(cluster["I"], 10);
	
	// Test HeVClusters
	
	for (int i = 30; i < 75; i++) {
		// Get the actual He and V amounts
		
		shared_ptr<Reactant> reactant = network->reactants->at(i);
		shared_ptr<HeVCluster> cluster = std::dynamic_pointer_cast<HeVCluster>(reactant);
		int actualHe = cluster->getSpeciesSize("He");
		int actualV = cluster->getSpeciesSize("V");
		
		// Get the amounts expected by Reactant::toClusterMap()
		
		std::map<std::string, int> species = network->toClusterMap(i);
		int expectedHe = species["He"];
		int expectedV = species["V"];
		
		BOOST_REQUIRE_EQUAL(actualHe, expectedHe);
		BOOST_REQUIRE_EQUAL(actualV, expectedV);
	}
}

BOOST_AUTO_TEST_CASE(checkReactants) {

	// Create the network
	auto psiNetwork = make_shared<PSIClusterReactionNetwork>();

	BOOST_FAIL("Unimplemented.");

	// Add a few He, V and I PSIClusters

	// Check the network

	// Add a couple of HeV and HeI clusters ("compounds")

	// Check the network

	// Try to add a couple of things that can't be there.

	// Make sure they didn't get added to the network.

	return;
}

BOOST_AUTO_TEST_CASE(checkProperties) {

	// Create the network
	auto psiNetwork = make_shared<PSIClusterReactionNetwork>();

	// Grab the map of properties from the network
	auto props = psiNetwork->getProperties();
	// Convert the property strings so we can use them
	auto numHeClusters = std::stoi(props["numHeClusters"]);
	auto numVClusters = std::stoi(props["numVClusters"]);
	auto numIClusters = std::stoi(props["numIClusters"]);
	auto numHeVClusters = std::stoi(props["numHeVClusters"]);
	auto numHeIClusters = std::stoi(props["numHeIClusters"]);
	auto maxHeVClusterSize = std::stoi(props["maxHeVClusterSize"]);
	auto maxHeIClusterSize = std::stoi(props["maxHeIClusterSize"]);
	auto maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	auto maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	auto maxIClusterSize = std::stoi(props["maxIClusterSize"]);

	// Check the properties
	BOOST_REQUIRE_EQUAL(0,numHeClusters);
	BOOST_REQUIRE_EQUAL(0,numVClusters);
	BOOST_REQUIRE_EQUAL(0,numIClusters);
	BOOST_REQUIRE_EQUAL(0,numHeVClusters);
	BOOST_REQUIRE_EQUAL(0,numHeIClusters);
	BOOST_REQUIRE_EQUAL(0,maxHeVClusterSize);
	BOOST_REQUIRE_EQUAL(0,maxHeIClusterSize);
	BOOST_REQUIRE_EQUAL(0,maxHeClusterSize);
	BOOST_REQUIRE_EQUAL(0,maxVClusterSize);
	BOOST_REQUIRE_EQUAL(0,maxIClusterSize);

	// Set a couple of properties
	psiNetwork->setProperty("rangePenalty","5");
	psiNetwork->setProperty("agility","d8");

	// Grab the properties afresh
	auto modifiedProps = psiNetwork->getProperties();

	// Check for the new properties
	auto rangePenalty = modifiedProps["rangePenalty"];
	auto agility = modifiedProps["agility"];
	BOOST_REQUIRE_EQUAL("5",rangePenalty);
	BOOST_REQUIRE_EQUAL("d8",agility);

	// Add a couple of clusters
	auto heCluster = make_shared<HeCluster>(5);
	psiNetwork->add(heCluster);
	auto heVCluster = make_shared<HeVCluster>(5,3);
	psiNetwork->add(heVCluster);

	// Grab the properties afresh
	auto propsWithClusters = psiNetwork->getProperties();
	numHeClusters = std::stoi(propsWithClusters["numHeClusters"]);
	maxHeClusterSize = std::stoi(propsWithClusters["maxHeClusterSize"]);
	numHeVClusters = std::stoi(propsWithClusters["numHeVClusters"]);
	maxHeVClusterSize = std::stoi(propsWithClusters["maxHeVClusterSize"]);

	// Check the properties again
	BOOST_REQUIRE_EQUAL(1,numHeClusters);
	BOOST_REQUIRE_EQUAL(1,numHeVClusters);
	BOOST_REQUIRE_EQUAL(5,maxHeClusterSize);
	BOOST_REQUIRE_EQUAL(8,maxHeVClusterSize);

	return;
}

BOOST_AUTO_TEST_CASE(checkNames) {

	// Create the network
	auto psiNetwork = make_shared<PSIClusterReactionNetwork>();

	BOOST_FAIL("Unimplemented.");

	// Create the names

	// Check the names of the regular cluster types

	// Check the names of the compound cluster types

	return;
}

BOOST_AUTO_TEST_SUITE_END()
