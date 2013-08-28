/*
 * SimpleReactionNetwork.cpp
 *
 *  Created on: Jun 29, 2013
 *      Author: bkj
 */

#include "SimpleReactionNetwork.h"
#include <PSICluster.h>
#include <HeCluster.h>
#include <VCluster.h>
#include <InterstitialCluster.h>
#include <HeVCluster.h>
#include <HeInterstitialCluster.h>
#include "SimpleReactionNetwork.h"
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using std::shared_ptr;
using namespace xolotlCore;
using namespace testUtils;

SimpleReactionNetwork::SimpleReactionNetwork() {

	// Hard code the size of the largest cluster
	int maxClusterSize = 10;
	int numClusters = maxClusterSize;

	// Add He clusters
	for (int numHe = 1; numHe <= maxClusterSize; numHe++) {
		// Create a He cluster with cluster size numHe
		std::shared_ptr<HeCluster> cluster(new HeCluster(numHe));

		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add vacancy clusters
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		// Create a He cluster with cluster size numV
		std::shared_ptr<VCluster> cluster(new VCluster(numV));

		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add interstitial clusters
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		// Create a He cluster with cluster size numI
		std::shared_ptr<InterstitialCluster> cluster(
				new InterstitialCluster(numI));

		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add HeV clusters, assuming that
	// numHe + numV <= maxMixedClusterSize

	int numHeVClusters = 0;

	for (int numV = 1; numV <= maxClusterSize; numV++) {
		for (int numHe = 1; numHe + numV <= maxClusterSize; numHe++) {
			// Create a HeVCluster with the current amount of He and V
			std::shared_ptr<HeVCluster> cluster(new HeVCluster(numHe, numV));

			reactants->push_back(cluster);
			numHeVClusters++;
		}
	}

	// Create the HeI clusters in this simple reaction network
	int numHeIClusters = 0;

	// Create all possible combinations of numHe and numI
	// clusters with numHe, numI < maxClusterSize
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		for (int numHe = 1; numHe + numI <= maxClusterSize; numHe++) {
			// Create the HeI cluster
			std::shared_ptr<HeInterstitialCluster> cluster(
					new HeInterstitialCluster(numHe, numI));

			// Add it to the reactants vector
			reactants->push_back(cluster);

			// Increment the number of
			// HeIClusters for the properties map
			numHeIClusters++;
		}
	}

	// Setup the properties map
	(*properties)["maxHeClusterSize"] = std::to_string(
			(long long) maxClusterSize);
	(*properties)["maxVClusterSize"] = std::to_string(
			(long long) maxClusterSize);
	(*properties)["maxIClusterSize"] = std::to_string(
			(long long) maxClusterSize);
	(*properties)["maxMixedClusterSize"] = std::to_string(
			(long long) maxClusterSize);

	(*properties)["numHeClusters"] = std::to_string((long long) numClusters);
	(*properties)["numVClusters"] = std::to_string((long long) numClusters);
	(*properties)["numIClusters"] = std::to_string((long long) numClusters);
	(*properties)["numHeVClusters"] = std::to_string(
			(long long) numHeVClusters);
	(*properties)["numHeIClusters"] = std::to_string(
			(long long) numHeIClusters);
}

SimpleReactionNetwork::~SimpleReactionNetwork() {
	// Nothing to do
}

/**
 * This operation creates a SimpleReactionNetwork and makes sure that it is
 * properly registered with the clusters it contains. This operation should
 * always be called instead of constructing a SimpleReactionNetwork manually.
 * @return The reaction network.
 */
std::shared_ptr<xolotlCore::ReactionNetwork> testUtils::getSimpleReactionNetwork() {
	// Create the network
	std::shared_ptr<xolotlCore::ReactionNetwork> network(
			new SimpleReactionNetwork());
	// Register the reaction network with its clusters
	for (int i = 0; i < network->reactants->size(); i++) {
		network->reactants->at(i)->setReactionNetwork(network);
	}

	// TEMPORARY
	// Print the reaction connectivity matrix

	for (auto reactantIt = network->reactants->begin();
			reactantIt != network->reactants->end(); reactantIt++) {
		std::shared_ptr<PSICluster> cluster = std::dynamic_pointer_cast
				< PSICluster > (*reactantIt);
		std::vector<int> conn = *cluster->getConnectivity();

		for (auto connIt = conn.begin(); connIt != conn.end(); connIt++) {
			printf("%s", *connIt ? "* " : "  ");
		}
		printf("\n");
	}

	return network;
}

/**
 * This operation returns a reactant with the given name and size if it
 * exists in the network or null if not.
 * @param name the name of the reactant
 * @param size the size of the reactant
 * @return A shared pointer to the reactant
 */
std::shared_ptr<Reactant> SimpleReactionNetwork::get(const std::string rName,
		const int size) {
	std::shared_ptr<Reactant> nullPtr;

	return nullPtr;
}

/**
 * This operation returns a compound reactant with the given name and size if it
 * exists in the network or null if not.
 * @param name the name of the compound reactant
 * @param sizes an array containing the sizes of each piece of the reactant
 * @return A shared pointer to the compound reactant
 */
std::shared_ptr<Reactant> SimpleReactionNetwork::getCompound(
		const std::string rName, const std::vector<int> sizes) {
	std::shared_ptr<Reactant> nullPtr;

	return nullPtr;
}

/**
 * This operation adds a reactant or a compound reactant to the network.
 * @param reactant The reactant that should be added to the network.
 */
void SimpleReactionNetwork::add(std::shared_ptr<Reactant> reactant) {
	return;
}

/**
 * This operation returns the names of the reactants in the network.
 * @return A vector with one each for each of the distinct reactant types
 * in the network.
 */
const std::vector<std::string> & SimpleReactionNetwork::getNames() {

	std::vector<std::string> vector;
	std::vector<std::string> & vectorRef = vector;

	return vector;
}

/**
 * This operation returns the names of the compound reactants in the
 * network.
 * @return A vector with one each for each of the distinct compound
 * reactant types in the network.
 */
const std::vector<std::string> & SimpleReactionNetwork::getCompoundNames() {
	std::vector<std::string> vector;
	std::vector<std::string> & vectorRef = vector;

	return vectorRef;
}

/**
 * This operation returns a map of the properties of this reaction network.
 * @return The map of properties that has been configured for this
 * ReactionNetwork.
 */
const std::map<std::string, std::string> & SimpleReactionNetwork::getProperties() {

	std::map<std::string, std::string> map;
	std::map<std::string, std::string> & mapRef = map;

	return mapRef;
}

/**
 * This operation sets a property with the given key to the specified value
 * for the network. ReactionNetworks may reserve the right to ignore this
 * operation for special key types.
 * @param key The key for the property
 * @param value The value to which the key should be set.
 */
void SimpleReactionNetwork::setProperty(const std::string key,
		const std::string value) {
	return;
}
