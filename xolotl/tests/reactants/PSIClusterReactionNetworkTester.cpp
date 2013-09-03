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
#include <VCluster.h>
#include <InterstitialCluster.h>
#include <HeInterstitialCluster.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the ReactionNetwork
 */BOOST_AUTO_TEST_SUITE(PSIReactionNetwork_testSuite)

BOOST_AUTO_TEST_CASE(checkCompositionCreation) {

	// Create the network
	auto psiNetwork = make_shared<PSIClusterReactionNetwork>();

	// Get an HeV cluster with sizes 5,10
	auto heVComp = psiNetwork->getCompositionVector(5, 10, 0);
	BOOST_REQUIRE_EQUAL(3, heVComp.size());
	BOOST_REQUIRE_EQUAL(5, heVComp[0]);
	BOOST_REQUIRE_EQUAL(10, heVComp[1]);
	BOOST_REQUIRE_EQUAL(0, heVComp[2]);

	// Get an interstitial cluster with size 99
	auto intComp = psiNetwork->getCompositionVector(0, 0, 99);
	BOOST_REQUIRE_EQUAL(3, intComp.size());
	BOOST_REQUIRE_EQUAL(0, intComp[0]);
	BOOST_REQUIRE_EQUAL(0, intComp[1]);
	BOOST_REQUIRE_EQUAL(99, intComp[2]);

	// Try to get something bad
	auto badComp = psiNetwork->getCompositionVector(-1, -8, -3);
	// Make sure that it gives you single species helium back
	BOOST_REQUIRE_EQUAL(3, badComp.size());
	BOOST_REQUIRE_EQUAL(1, badComp[0]);
	BOOST_REQUIRE_EQUAL(0, badComp[1]);
	BOOST_REQUIRE_EQUAL(0, badComp[2]);

}

BOOST_AUTO_TEST_CASE(checkReactants) {

	// Create the network
	auto psiNetwork = make_shared<PSIClusterReactionNetwork>();

	// Add a few He, V and I PSIClusters
	auto heCluster = make_shared<HeCluster>(10);
	auto vCluster = make_shared<VCluster>(4);
	auto interstitialCluster = make_shared<InterstitialCluster>(48);
	psiNetwork->add(heCluster);
	psiNetwork->add(vCluster);
	psiNetwork->add(interstitialCluster);

	// Check the network, He first
	auto retHeCluster = std::dynamic_pointer_cast<PSICluster>(
			psiNetwork->get("He", 10));
	BOOST_REQUIRE(retHeCluster);
	BOOST_REQUIRE_EQUAL("He", retHeCluster->getName());
	BOOST_REQUIRE_EQUAL(10, retHeCluster->getSize());
	// V
	auto retVCluster = std::dynamic_pointer_cast<PSICluster>(
			psiNetwork->get("V", 4));
	BOOST_REQUIRE(retVCluster);
	BOOST_REQUIRE_EQUAL(4, retVCluster->getSize());
	BOOST_REQUIRE_EQUAL("V", retVCluster->getName());
	// I
	auto retICluster = std::dynamic_pointer_cast<PSICluster>(
			psiNetwork->get("I", 48));
	BOOST_REQUIRE(retICluster);
	BOOST_REQUIRE_EQUAL(48, retICluster->getSize());
	BOOST_REQUIRE_EQUAL("I", retICluster->getName());

	// Add a couple of HeV and HeI clusters ("compounds")
	auto heVCluster = make_shared<HeVCluster>(5, 3);
	auto heICluster = make_shared<HeInterstitialCluster>(8, 8);
	psiNetwork->add(heVCluster);
	psiNetwork->add(heICluster);

	// Check the network, start with HeV
	std::vector<int> sizes = psiNetwork->getCompositionVector(5, 3, 0);
	std::shared_ptr<PSICluster> retHeVCluster = std::dynamic_pointer_cast<
			PSICluster>(psiNetwork->getCompound("HeV", sizes));
	BOOST_REQUIRE_EQUAL(8, retHeVCluster->getSize());
	BOOST_REQUIRE_EQUAL("HeV", retHeVCluster->getName());
	// Now do HeI
	sizes = psiNetwork->getCompositionVector(8, 0, 8);
	auto retHeICluster = std::dynamic_pointer_cast<PSICluster>(
			psiNetwork->getCompound("HeI", sizes));
	BOOST_REQUIRE_EQUAL(16, retHeICluster->getSize());
	BOOST_REQUIRE_EQUAL("HeI", retHeICluster->getName());

	// Check the getter for all reactants
	auto clusters = psiNetwork->getAll();
	BOOST_REQUIRE_EQUAL(5, clusters->size());

	// Check the size of the network
	BOOST_REQUIRE_EQUAL(5, psiNetwork->size());

	// Check the cluster ids. All we can do is check that they are between 1
	// and 5. Start with the He cluster
	int id = psiNetwork->getReactantId(*retHeCluster);
	BOOST_REQUIRE(id > 0 && id <= 5);
	// V
	id = psiNetwork->getReactantId(*retVCluster);
	BOOST_REQUIRE(id > 0 && id <= 5);
	// I
	id = psiNetwork->getReactantId(*retICluster);
	BOOST_REQUIRE(id > 0 && id <= 5);
	// HeV
	id = psiNetwork->getReactantId(*retHeVCluster);
	BOOST_REQUIRE(id > 0 && id <= 5);
	// HeI
	id = psiNetwork->getReactantId(*retHeICluster);
	BOOST_REQUIRE(id > 0 && id <= 5);

	// Add a whole bunch of HeV clusters to make sure that the network can
	// handle large numbers of them properly.
	int counter = 0;
	int maxClusterSize = 10;
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		for (int numHe = 1; numHe + numV <= maxClusterSize; numHe++) {
			std::shared_ptr<HeVCluster> cluster(new HeVCluster(numHe, numV));
			psiNetwork->add(cluster);
			counter++;
		}
	}
	BOOST_TEST_MESSAGE("Added " << counter << " HeV clusters");

	// Add a whole bunch of HeI clusters to make sure that the network can
	// handle large numbers of them properly too. Use a different max cluster
	// size to throw a little variability at the network.
	counter = 0;
	maxClusterSize = 9;
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		for (int numHe = 1; numHe + numI <= maxClusterSize; numHe++) {
			std::shared_ptr<HeInterstitialCluster> cluster(
					new HeInterstitialCluster(numHe, numI));
			psiNetwork->add(cluster);
			counter++;
		}
	}
	BOOST_TEST_MESSAGE("Added " << counter << " HeI clusters");

	// Make sure that everything was added
	auto reactants = psiNetwork->getAll();
	BOOST_REQUIRE_EQUAL(85,reactants->size());
	// Get the clusters by type and check them. Start with He.
	reactants = psiNetwork->getAll("He");
	BOOST_REQUIRE_EQUAL(1, reactants->size());
	BOOST_REQUIRE_EQUAL("He", reactants->at(0)->getName());
	// V
	reactants = psiNetwork->getAll("V");
	BOOST_REQUIRE_EQUAL(1, reactants->size());
	BOOST_REQUIRE_EQUAL("V", reactants->at(0)->getName());
	// I
	reactants = psiNetwork->getAll("I");
	BOOST_REQUIRE_EQUAL(1, reactants->size());
	BOOST_REQUIRE_EQUAL("I", reactants->at(0)->getName());
	// HeV
	reactants = psiNetwork->getAll("HeV");
	BOOST_REQUIRE_EQUAL(45, reactants->size());
	BOOST_REQUIRE_EQUAL("HeV", reactants->at(0)->getName());
	// HeI
	reactants = psiNetwork->getAll("HeI");
	BOOST_REQUIRE_EQUAL(37, reactants->size());
	BOOST_REQUIRE_EQUAL("HeI", reactants->at(0)->getName());

	// Try to get something that obviously isn't there

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
	BOOST_REQUIRE_EQUAL(0, numHeClusters);
	BOOST_REQUIRE_EQUAL(0, numVClusters);
	BOOST_REQUIRE_EQUAL(0, numIClusters);
	BOOST_REQUIRE_EQUAL(0, numHeVClusters);
	BOOST_REQUIRE_EQUAL(0, numHeIClusters);
	BOOST_REQUIRE_EQUAL(0, maxHeVClusterSize);
	BOOST_REQUIRE_EQUAL(0, maxHeIClusterSize);
	BOOST_REQUIRE_EQUAL(0, maxHeClusterSize);
	BOOST_REQUIRE_EQUAL(0, maxVClusterSize);
	BOOST_REQUIRE_EQUAL(0, maxIClusterSize);

	// Set a couple of properties
	psiNetwork->setProperty("rangePenalty", "5");
	psiNetwork->setProperty("agility", "d8");

	// Grab the properties afresh
	auto modifiedProps = psiNetwork->getProperties();

	// Check for the new properties
	auto rangePenalty = modifiedProps["rangePenalty"];
	auto agility = modifiedProps["agility"];
	BOOST_REQUIRE_EQUAL("5", rangePenalty);
	BOOST_REQUIRE_EQUAL("d8", agility);

	// Add a couple of clusters
	auto heCluster = make_shared<HeCluster>(5);
	psiNetwork->add(heCluster);
	auto heVCluster = make_shared<HeVCluster>(5, 3);
	psiNetwork->add(heVCluster);

	// Grab the properties afresh
	auto propsWithClusters = psiNetwork->getProperties();
	numHeClusters = std::stoi(propsWithClusters["numHeClusters"]);
	maxHeClusterSize = std::stoi(propsWithClusters["maxHeClusterSize"]);
	numHeVClusters = std::stoi(propsWithClusters["numHeVClusters"]);
	maxHeVClusterSize = std::stoi(propsWithClusters["maxHeVClusterSize"]);

	// Check the properties again
	BOOST_REQUIRE_EQUAL(1, numHeClusters);
	BOOST_REQUIRE_EQUAL(1, numHeVClusters);
	BOOST_REQUIRE_EQUAL(5, maxHeClusterSize);
	BOOST_REQUIRE_EQUAL(8, maxHeVClusterSize);

	return;
}

BOOST_AUTO_TEST_CASE(checkNames) {

	// Create the network
	auto psiNetwork = make_shared<PSIClusterReactionNetwork>();

	// Check the names of the regular cluster types. Use a simple counting
	// system to look over the list since there is no way to check exact
	// containment with a vector.
	auto names = psiNetwork->getNames();
	int marker = 0;
	for (int i = 0; i < names.size(); i++) {
		if (names[i] == "He")
			++marker;
		else if (names[i] == "V")
			++marker;
		else if (names[i] == "I")
			++marker;
	}
	BOOST_REQUIRE_EQUAL(3, marker);
	BOOST_REQUIRE_EQUAL(marker, names.size());

	// Check the names of the compound cluster types. Use the same counting
	// system as above.
	auto compoundNames = psiNetwork->getCompoundNames();
	marker = 0;
	for (int i = 0; i < compoundNames.size(); i++) {
		if (compoundNames[i] == "HeV")
			++marker;
		else if (compoundNames[i] == "HeI")
			++marker;
	}
	BOOST_REQUIRE_EQUAL(2, marker);
	BOOST_REQUIRE_EQUAL(marker, compoundNames.size());

	return;
}

/**
 * This operation tests the copy constructor.
 */BOOST_AUTO_TEST_CASE(checkCopying) {

	PSIClusterReactionNetwork network;

	// Add a reactant
	std::shared_ptr<Reactant> heCluster(new HeCluster(1));
	heCluster->setConcentration(50.0);
	network.add(heCluster);

	// Copy the network
	PSIClusterReactionNetwork networkCopy = network;

	// Check that the ReactionNetwork fields are copied
	auto properties = network.getProperties();
	auto copiedProperties = networkCopy.getProperties();
	BOOST_REQUIRE_EQUAL(properties.size(), copiedProperties.size());
	BOOST_REQUIRE_EQUAL(properties["numHeClusters"],
			copiedProperties["numHeClusters"]);

	// Check that changing the concentration of a copy does not update the
	// original. Start by updating the copy.
	auto copiedHeCluster = networkCopy.get("He", 1);
	copiedHeCluster->setConcentration(7.0);
	BOOST_REQUIRE_CLOSE(7.0, copiedHeCluster->getConcentration(), 1.0e-5);
	// Make sure the original wasn't changed.
	BOOST_REQUIRE_CLOSE(50.0, heCluster->getConcentration(), 1.0e-5);

	// Check the size of the network
	BOOST_REQUIRE_EQUAL(1, networkCopy.size());

	return;
}
BOOST_AUTO_TEST_SUITE_END()
