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
#include "SimpleReactionNetwork.h"
#include <HeCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using std::shared_ptr;
using namespace xolotlCore;
using namespace testUtils;

/**
 * This suite is responsible for testing the HeCluster.
 */
BOOST_AUTO_TEST_SUITE(HeCluster_testSuite)

/**
 * This operation checks the ability of the HeCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	std::vector<shared_ptr<Reactant>> &reactants = *network->reactants;
	std::map<std::string, std::string> &props = *network->properties;
	
	// Store network properties into local variables
	
	int numHeClusters = std::stoi(props["numHeClusters"]);
	int numVClusters = std::stoi(props["numVClusters"]);
	int numIClusters = std::stoi(props["numIClusters"]);
	int numHeVClusters = std::stoi(props["numHeVClusters"]);
	int maxMixedClusterSize = std::stoi(props["maxMixedClusterSize"]);
	
	// Offsets determine the position of the first type of cluster in the
	// reactants array
	
	int offsetV = numHeClusters;
	int offsetI = offsetV + numVClusters;
	int offsetHeV = offsetI + numIClusters;
	
	// Check the connectivity of the first He reactant (numHe=1)
	
	{
		std::vector<int> connectivityArray = reactants.at(0)->getConnectivity();
		int numHe = 1;
		
		// Check the connectivity for He, V, and I
		
		int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
			
			// V
			1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
			
			// I
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0
		};
		
		for (int i = 0; i < 30; i++) {
			BOOST_REQUIRE_EQUAL(connectivityArray.at(i), connectivityExpected[i]);
		}
		
		// Check the connectivity for HeV
		
		for (int i = 30; i < reactants.size(); i++) {
			// Build a species map
			std::map<std::string, int> speciesMap = network->toClusterMap(i);
			int numHeOther = speciesMap["He"];
			int numVOther = speciesMap["V"];
			
			bool connected = numHe + numHeOther + numVOther <= maxMixedClusterSize;
			BOOST_REQUIRE_EQUAL(connectivityArray.at(i), (int) connected);
		}
	}
	
	
	// Check the connectivity of the 6th He reactant (numHe=6)
	
	{
		std::vector<int> connectivityArray = reactants.at(5)->getConnectivity();
		int numHe = 6;
		
		// Check the connectivity for He, V, and I
		
		int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			
			// V
			1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			
			// I
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0
		};
		
		for (int i = 0; i < 30; i++) {
			BOOST_REQUIRE_EQUAL(connectivityArray.at(i), connectivityExpected[i]);
		}
		
		// Check the connectivity for HeV
		
		for (int i = 30; i < network->reactants->size(); i++) {
			// Build a species map
			std::map<std::string, int> speciesMap = network->toClusterMap(i);
			int numHeOther = speciesMap["He"];
			int numVOther = speciesMap["V"];
			
			bool connected = numHe + numHeOther + numVOther <= maxMixedClusterSize;
			BOOST_REQUIRE_EQUAL(connectivityArray.at(i), (int) connected);
		}
	}
}

/**
 * This operation checks the HeCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
}

/**
 * This operation checks the reaction radius for HeCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {

	std::vector<std::shared_ptr<HeCluster>> clusters;
	std::shared_ptr<HeCluster> cluster;

	double expectedRadii[] = { 0.3, 0.3748419767, 0.4273418681, 0.4691369586,
			0.5044313198, 0.5352826768, 0.5628704922, 0.5879411911,
			0.6110006225, 0.6324092998 };

	for (int i = 1; i <= 10; i++) {
		cluster = std::shared_ptr<HeCluster>(new HeCluster(i));
		BOOST_CHECK_CLOSE(expectedRadii[i-1], cluster->getReactionRadius(), .000001);
	}
}


BOOST_AUTO_TEST_SUITE_END()

