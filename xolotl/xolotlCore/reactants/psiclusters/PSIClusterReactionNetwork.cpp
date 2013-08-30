#include "PSIClusterReactionNetwork.h"
#include "PSICluster.h"
#include <iostream>
#include <algorithm>

using namespace xolotlCore;
using std::shared_ptr;

void PSIClusterReactionNetwork::setDefaultPropsAndNames() {

	// Initialize default properties
	(*properties)["reactionsEnabled"] = "true";
	(*properties)["dissociationsEnabled"] = "true";
	(*properties)["numHeClusters"] = "0";
	(*properties)["numVClusters"] = "0";
	(*properties)["numIClusters"] = "0";
	(*properties)["numHeVClusters"] = "0";
	(*properties)["numHeIClusters"] = "0";
	(*properties)["maxHeClusterSize"] = "0";
	(*properties)["maxVClusterSize"] = "0";
	(*properties)["maxIClusterSize"] = "0";
	(*properties)["maxHeVClusterSize"] = "0";
	(*properties)["maxHeIClusterSize"] = "0";

	// Set the reactant names
	names.push_back("He");
	names.push_back("V");
	names.push_back("I");
	// Set the compound reactant names
	compoundNames.push_back("HeV");
	compoundNames.push_back("HeI");

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork() :
		ReactionNetwork() {

	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork(
		const PSIClusterReactionNetwork &other) :
		ReactionNetwork(other) {

	// Local Declarations
	int numNames, i;

	// Copy the names
	names = other.getNames();
	// Copy the compound names
	compoundNames = other.getCompoundNames();
	// Reset the properties table so that it can be properly updated when the
	// network is filled.
	setDefaultPropsAndNames();
	// Get all of the reactants from the other network and add them to this one
	auto reactants = other.getAll();
	for (i = 0; i < reactants->size(); i++) {
		add(reactants->at(i)->clone());
	}

}

std::map<std::string, int> PSIClusterReactionNetwork::toClusterMap(
		int index) const {

	// This functionality is implemented in PSICluster.
	// Instead of redefining the map for all indices, we simply obtain
	// the cluster map from the reactant in the reactants vector at the
	// position of the index.
	shared_ptr<PSICluster> cluster = std::dynamic_pointer_cast<PSICluster>(
			reactants->at(index));
	return cluster->getClusterMap();
}

int PSIClusterReactionNetwork::toClusterIndex(
		std::map<std::string, int> clusterMap) const {

	// If the clusterMap doesn't have one of the following keys,
	// the value will be set to zero.
	int numHe = 0;
	int numV = 0;
	int numI = 0;
	int indexOffset = 0;
	int numReactants = reactants->size();
	int finalIndex = -1;
	std::map<std::string, int> otherClusterMap;
	std::shared_ptr<PSICluster> cluster;
	// Convert the property strings so we can use them
	int numHeClusters = std::stoi((*properties)["numHeClusters"]);
	int numVClusters = std::stoi((*properties)["numVClusters"]);
	int numIClusters = std::stoi((*properties)["numIClusters"]);
	int numHeVClusters = std::stoi((*properties)["numHeVClusters"]);
	int numHeIClusters = std::stoi((*properties)["numHeIClusters"]);
	int maxMixedClusterSize = std::stoi((*properties)["maxMixedClusterSize"]);
	int numSpecies = 0;

	// Update the indices
	numHe = std::max(clusterMap["He"], 0);
	numV = std::max(clusterMap["V"], 0);
	numI = std::max(clusterMap["I"], 0);
	numSpecies = (numHe > 0) + (numV > 0) + (numI > 0);

	// Cluster should either contain one or two types of species.
	if (numSpecies == 0) {
		throw std::string("Cluster map contains no species");
	} else if (numSpecies == 1) {
		// Single species
		if (numHe) {
			finalIndex = numHe - 1;
		} else if (numV) {
			finalIndex = numV + numHeClusters - 1;
		} else if (numI) {
			finalIndex = numI + numHeClusters + numVClusters - 1;
		}
	} else if (numSpecies == 2) {
		// HeVCluster
		if (numHe && numV) {
			indexOffset = numHeClusters + numVClusters + numIClusters;
			// Do a linear search across the set of HeV clusters to find the
			// cluster. FIXME - THIS IS SLOW FOR LARGE N!
//			for (int i = indexOffset; i < indexOffset + numHeVClusters; i++) {
//				cluster = std::dynamic_pointer_cast<PSICluster>(reactants->at(i));
//				otherClusterMap = cluster->getClusterMap();
//				if (otherClusterMap.at("He") == numHe && otherClusterMap.at("V") == numV) {
//					finalIndex = i;
//					break;
//				}
//			}
			// FIXME FIXME FIXME!!!!!!!!
			finalIndex = 1;
		}
		// ----- HeICluster -----
		// Increment the offset by the number of HeVClusters
		// indexOffset += maxMixedClusterSize * (maxMixedClusterSize - 1) / 2;
		else if (numHe && numI && numHeIClusters > 0) {
			indexOffset = numHeClusters + numVClusters + numIClusters
					+ numHeVClusters;
			// Closed form for converting a top-left triangle grid
			// to an index
			int index = (numI - 1) * maxMixedClusterSize - numI * (numI - 1) / 2
					+ numHe - 1;
			finalIndex = indexOffset + index;
		}
	}

	return finalIndex;
}

/**
 * This operation returns a reactant with the given name and size if it
 * exists in the network or null if not.
 * @param name the name of the reactant
 * @param size the size of the reactant
 * @return A shared pointer to the reactant
 */
std::shared_ptr<Reactant> PSIClusterReactionNetwork::get(
		const std::string rName, const int size) const {

	// Local Declarations
	std::map<std::string, int> composition;
	std::shared_ptr<PSICluster> retReactant = std::shared_ptr<PSICluster>();

	// Setup the composition map to default values
	composition["He"] = 0;
	composition["V"] = 0;
	composition["I"] = 0;

	// Only pull the reactant if the name is valid
	if (rName == "He" || rName == "V" || rName == "I") {
		composition[rName] = size;
		retReactant = singleSpeciesMap.at(composition);
	}

	return retReactant;
}

/**
 * This operation returns a compound reactant with the given name and size if it
 * exists in the network or null if not.
 * @param rName the name of the compound reactant
 * @param sizes an array containing the sizes of each piece of the reactant
 * @return A shared pointer to the compound reactant
 */
std::shared_ptr<Reactant> PSIClusterReactionNetwork::getCompound(
		const std::string rName, const std::vector<int> sizes) const {
	std::shared_ptr<Reactant> nullPtr;

	// Local Declarations
	std::map<std::string, int> composition;
	std::shared_ptr<PSICluster> retReactant;

	// Setup the composition map to default values
	composition["He"] = 0;
	composition["V"] = 0;
	composition["I"] = 0;

	// Only pull the reactant if the name is valid and there are enough sizes
	// to fill the composition.
	if ((rName == "HeV" || rName == "HeI") && sizes.size()) {
		composition["He"] = sizes[0];
		composition["V"] = sizes[1];
		composition["I"] = sizes[2];
		retReactant = mixedSpeciesMap.at(composition);
	}

	return retReactant;
}

/**
 * This operation returns all reactants in the network without regard for
 * their composition or whether they are compound reactants. The list may
 * or may not be ordered and the decision is left to implementers.
 * @return The list of all of the reactants in the network
 */
std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>PSIClusterReactionNetwork::getAll() const {

	// Local Declarations
	std::shared_ptr<std::vector<std::shared_ptr<Reactant>>> allReactants(new std::vector<std::shared_ptr<Reactant> >);

	// Load the single-species clusters
	for (auto it=singleSpeciesMap.begin(); it!=singleSpeciesMap.end(); ++it) {
		allReactants->push_back(it->second);
	}

	// Load the mixed-species clusters
	for (auto it=mixedSpeciesMap.begin(); it!=mixedSpeciesMap.end(); ++it) {
		allReactants->push_back(it->second);
	}

	return allReactants;
}

/**
 * This operation adds a reactant or a compound reactant to the network.
 * @param reactant The reactant that should be added to the network.
 */
void PSIClusterReactionNetwork::add(std::shared_ptr<Reactant> reactant) {

	// Local Declarations
	int numHe = 0, numV = 0, numI = 0;
	bool isMixed = false;
	std::string numClusterKey, clusterSizeKey;

	// Only add a complete reactant
	if (reactant != NULL) {

		// Get the composition
		auto composition = reactant->getComposition();

		// Get the cluster sizes
		numHe = composition.at("He");
		numV = composition.at("V");
		numI = composition.at("I");

		// Determine if the cluster is a compound. If there is more than one
		// type, then the check below will sum to greater than one and we know
		// that we have a mixed cluster.
		isMixed = ((numHe > 0) + (numV > 0) + (numI > 0)) > 1;

		// Add the compound or regular reactant.
		if (isMixed) {
			mixedSpeciesMap[composition] =
					std::dynamic_pointer_cast<PSICluster>(reactant);
			// Figure out whether we have HeV or HeI and set the keys
			if (numV > 0) {
				numClusterKey = "numHeVClusters";
				clusterSizeKey = "maxHeVClusterSize";
			} else {
				numClusterKey = "numHeIClusters";
				clusterSizeKey = "maxHeIClusterSize";
			}
		} else if (!isMixed) {
			singleSpeciesMap[composition] =
					std::dynamic_pointer_cast<PSICluster>(reactant);
			// Figure out whether we have He, V or I and set the keys
			if (numHe > 0) {
				numClusterKey = "numHeClusters";
				clusterSizeKey = "maxHeClusterSize";
			} else if (numV > 0) {
				numClusterKey = "numHeVClusters";
				clusterSizeKey = "maxVClusterSize";
			} else {
				numClusterKey = "numIClusters";
				clusterSizeKey = "maxIClusterSize";
			}
		}
		// Increment the number of total clusters of this type
		int numClusters = std::stoi(properties->at(numClusterKey));
		numClusters++;
		(*properties)[numClusterKey] = std::to_string((long long) numClusters);
		// Increment the max cluster size key
		int maxSize = std::stoi(properties->at(clusterSizeKey));
		int clusterSize = numHe + numV + numI;
		maxSize = std::max(clusterSize, maxSize);
		(*properties)[clusterSizeKey] = std::to_string((long long) maxSize);

	}

	return;
}

/**
 * This operation returns the names of the reactants in the network.
 * @return A vector with one each for each of the distinct reactant types
 * in the network.
 */
const std::vector<std::string> & PSIClusterReactionNetwork::getNames() const {

	return names;
}

/**
 * This operation returns the names of the compound reactants in the
 * network.
 * @return A vector with one each for each of the distinct compound
 * reactant types in the network.
 */
const std::vector<std::string> & PSIClusterReactionNetwork::getCompoundNames() const {

	return compoundNames;
}

/**
 * This operation returns a map of the properties of this reaction network.
 * @return The map of properties that has been configured for this
 * ReactionNetwork.
 */
const std::map<std::string, std::string> & PSIClusterReactionNetwork::getProperties() {

	return *properties;
}

/**
 * This operation sets a property with the given key to the specified value
 * for the network. ReactionNetworks may reserve the right to ignore this
 * operation for special key types, most especially those that they manage on
 * their own.
 *
 * Calling this operation will ignore all of the published properties that are
 * configured by default.
 * @param key The key for the property
 * @param value The value to which the key should be set.
 */
void PSIClusterReactionNetwork::setProperty(const std::string key,
		const std::string value) {

	// Check the keys and value before trying to set the property
	if (!key.empty() && !value.empty() && key != "numHeClusters"
			&& key != "numVClusters" && key != "numIClusters"
			&& key != "maxHeClusterSize" && key != "maxVClusterSize"
			&& key != "maxIClusterSize" && key != "maxHeVClusterSize"
			&& key != "maxHeIClusterSize") {
		// Add the property if it made it through that!
		(*properties)[key] = value;
	}

	return;
}
