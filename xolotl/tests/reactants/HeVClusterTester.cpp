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
#include <HeVCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;


/**
 * This suite is responsible for testing the HeVCluster.
 */
BOOST_AUTO_TEST_SUITE(HeVCluster_testSuite)


BOOST_AUTO_TEST_CASE(getSpeciesSize) {
	HeVCluster cluster(4, 5);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("He"), 4);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("V"), 5);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("I"), 0);
}

/**
 * This operation checks the ability of the HeVCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	int maxClusterSize = 10;
	int numClusters = maxClusterSize;
	
	shared_ptr<vector<shared_ptr<Reactant>>> reactants = network->reactants;
	shared_ptr<map<string, string>> props = network->properties;

	// Write the cluster information to stdout
	BOOST_TEST_MESSAGE("Sizes of clusters in network:");
	
	for (auto reactantIt = reactants->begin(); reactantIt != reactants->end(); reactantIt++) {
		// Write the size of the psi cluster to stdout
		shared_ptr<PSICluster> psiCluster = static_pointer_cast<PSICluster>(*reactantIt);
		BOOST_TEST_MESSAGE(psiCluster->getSize());
	}
	
	BOOST_TEST_MESSAGE("Maximum He Cluster Size = " << (*props)["maxHeClusterSize"]);
	BOOST_TEST_MESSAGE("Maximum V Cluster Size = " << (*props)["maxVClusterSize"]);
	BOOST_TEST_MESSAGE("Number of He clusters = " << (*props)["numHeClusters"]);
	BOOST_TEST_MESSAGE("Number of V clusters = " << (*props)["numVClusters"]);
	BOOST_TEST_MESSAGE("Number of mixed clusters = " << (*props)["numMixedClusters"]);

	// Get the connectivity of the 20th HeV cluster (index 59).
	shared_ptr<Reactant> reactant = reactants->at(3 * numClusters + 20 - 1);
	shared_ptr<HeVCluster> cluster = dynamic_pointer_cast<HeVCluster>(reactant);
	
	vector<int> connectivityArray = cluster->getConnectivity();
	
	// The connectivity array should be the same size as the reactants array
	BOOST_TEST_MESSAGE("Connectivity Array Size = " << connectivityArray.size());
	BOOST_REQUIRE_EQUAL(connectivityArray.size(), reactants->size());
	
	// The HeVCluster should have 6 He and 3 V, by the ordering of
	// the SimpleReactionNetwork.
	
	BOOST_REQUIRE_EQUAL(cluster->getSpeciesSize("He"), 6);
	BOOST_REQUIRE_EQUAL(cluster->getSpeciesSize("V"), 3);
	
	// The following are reactions of the HeVCluster
	
	// HeV[x, y] + He[z] --> HeV[x + z, y]
	// HeV[x, y] + V[1]  --> HeV[x, y + 1]
	// HeV[x, y] + I[z]  --> HeV[x, y - z]
	
	int connectivityExpected[] = {
		// He
		1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
		
		// V
		1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		
		// I
		1, 1, 0, 0, 0, 0, 0, 0, 0, 0
	};
	
	for (int i = 0; i < 30; i++) {
		BOOST_REQUIRE_EQUAL(connectivityArray.at(i), connectivityExpected[i]);
	}
	
	// Everything else should be 0
	int connectivitySize = connectivityArray.size();
	for (int i = 30; i < connectivitySize; i++) {
		BOOST_REQUIRE_EQUAL(connectivityArray.at(i), 0);
	}
}

BOOST_AUTO_TEST_CASE(checkGetFlux) {
	/*shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	int maxClusterSize = 10;
	int numClusters = maxClusterSize;

	shared_ptr<vector<shared_ptr<Reactant>>> reactants = network->reactants;
	shared_ptr<map<string, string>> props = network->properties;

	for (int i = 0; i < reactants->size(); i++) {
		reactants->at(i)->setConcentration(2.0*i);
		(dynamic_pointer_cast<PSICluster>(reactants->at(i)))->setMigrationEnergy(1.0);
		(dynamic_pointer_cast<PSICluster>(reactants->at(i)))->setDiffusionFactor(1.0);
	}
	// Get the connectivity of the 20th HeV cluster (index 59).
	shared_ptr<Reactant> reactant = reactants->at(3 * numClusters + 20 - 1);
	shared_ptr<HeVCluster> cluster = dynamic_pointer_cast<HeVCluster>(reactant);

	double flux = cluster->getProductionFlux(1.0);

	std::cout << "Flux is " << flux << "\n";*/
}

BOOST_AUTO_TEST_SUITE_END()
