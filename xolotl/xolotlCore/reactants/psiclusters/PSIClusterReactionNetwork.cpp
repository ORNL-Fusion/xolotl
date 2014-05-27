#include "PSIClusterReactionNetwork.h"
#include "PSICluster.h"
#include <HandlerRegistryFactory.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <Constants.h>

using namespace xolotlCore;
using std::shared_ptr;

void PSIClusterReactionNetwork::setDefaultPropsAndNames() {

	// Shared pointers for the cluster type map
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>>heVector =
			std::make_shared<std::vector<std::shared_ptr<Reactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>> vVector
		= std::make_shared<std::vector<std::shared_ptr<Reactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>> iVector
		= std::make_shared<std::vector<std::shared_ptr<Reactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>> heVVector
		= std::make_shared<std::vector<std::shared_ptr<Reactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>> heIVector
		= std::make_shared<std::vector<std::shared_ptr<Reactant>>>();

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

	// Initialize the current and last size to 0
	networkSize = 0;
	// Set the reactant names
	names.push_back(heType);
	names.push_back(vType);
	names.push_back(iType);
	// Set the compound reactant names
	compoundNames.push_back("HeV");
	compoundNames.push_back("HeI");

	// Setup the cluster type map
	clusterTypeMap[heType] = heVector;
	clusterTypeMap[vType] = vVector;
	clusterTypeMap[iType] = iVector;
	clusterTypeMap["HeV"] = heVVector;
	clusterTypeMap["HeI"] = heIVector;

	// Setup the vector to hold all of the reactants
	allReactants = make_shared<std::vector<Reactant *>>();

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork() :
		ReactionNetwork() {

	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork(registry) {

	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork(
		const PSIClusterReactionNetwork &other) :
		ReactionNetwork(other) {

	// The size and ids do not need to be copied. They will be fixed when the
	// reactants are added.

	// Copy the names
	names = other.getNames();
	// Copy the compound names
	compoundNames = other.getCompoundNames();
	// Reset the properties table so that it can be properly updated when the
	// network is filled.
	setDefaultPropsAndNames();
	// Get all of the reactants from the other network and add them to this one
	// Load the single-species clusters. Calling getAll() will not work because
	// it is not const.
	std::vector<std::shared_ptr<Reactant> > reactants;
	for (auto it = other.singleSpeciesMap.begin();
			it != other.singleSpeciesMap.end(); ++it) {
		reactants.push_back(it->second);
	}
	// Load the mixed-species clusters
	for (auto it = other.mixedSpeciesMap.begin();
			it != other.mixedSpeciesMap.end(); ++it) {
		reactants.push_back(it->second);
	}
	for (int i = 0; i < reactants.size(); i++) {
		add(reactants[i]->clone());
	}

}

/**
 * This operation sets the temperature at which the reactants currently
 * exists. It calls setTemperature() on each reactant.
 *
 * This is the simplest way to set the temperature for all reactants is to
 * call the ReactionNetwork::setTemperature() operation.
 *
 * @param temp The new temperature
 */
void PSIClusterReactionNetwork::setTemperature(double temp) {
	// Set the temperature
	temperature = temp;

	// Update the temperature for all of the clusters
	int networkSize = size();
	for (int i = 0; i < networkSize; i++) {
		allReactants->at(i)->setTemperature(temp);
	}

	return;
}

/**
 * This operation returns the temperature at which the cluster currently exists.
 * @return The temperature.
 */
double PSIClusterReactionNetwork::getTemperature() const {
	return temperature;
}

/**
 * This operation returns a reactant with the given name and size if it
 * exists in the network or null if not.
 * @param type the name of the reactant
 * @param size the size of the reactant
 * @return A shared pointer to the reactant
 */
Reactant * PSIClusterReactionNetwork::get(const std::string type,
		const int size) const {

	// Local Declarations
	static std::map<std::string, int> composition = { { heType, 0 },
			{ vType, 0 }, { iType, 0 } };
	std::shared_ptr<PSICluster> retReactant;

	// Setup the composition map to default values
	composition[heType] = 0;
	composition[vType] = 0;
	composition[iType] = 0;

	// Only pull the reactant if the name and size are valid
	if ((type == heType || type == vType || type == iType) && size >= 1) {
		composition[type] = size;
		//std::string encodedName = PSICluster::encodeCompositionAsName(composition);
		// Make sure the reactant is in the map
		if (singleSpeciesMap.count(composition)) {
			retReactant = singleSpeciesMap.at(composition);
		}
	}

	return (Reactant *) retReactant.get();
}

/**
 * This operation returns a compound reactant with the given name and size if it
 * exists in the network or null if not.
 * @param type the type of the compound reactant
 * @param sizes an array containing the sizes of each piece of the reactant
 * @return A shared pointer to the compound reactant
 */
Reactant * PSIClusterReactionNetwork::getCompound(const std::string type,
		const std::vector<int> sizes) const {

	// Local Declarations
	static std::map<std::string, int> composition = { { heType, 0 },
			{ vType, 0 }, { iType, 0 } };
	std::shared_ptr<PSICluster> retReactant;

	// Setup the composition map to default values
	composition[heType] = 0;
	composition[vType] = 0;
	composition[iType] = 0;

	// Only pull the reactant if the name is valid and there are enough sizes
	// to fill the composition.
	if ((type == "HeV" || type == "HeI") && sizes.size() == 3) {
		composition[heType] = sizes[0];
		composition[vType] = sizes[1];
		composition[iType] = sizes[2];
		// Make sure the reactant is in the map
		if (mixedSpeciesMap.count(composition)) {
			retReactant = mixedSpeciesMap.at(composition);
		}
	}

	return (Reactant *) retReactant.get();
}

/**
 * This operation returns all reactants in the network without regard for
 * their composition or whether they are compound reactants. The list may
 * or may not be ordered and the decision is left to implementers.
 * @return The list of all of the reactants in the network
 */
const std::shared_ptr<std::vector<Reactant *>> & PSIClusterReactionNetwork::getAll() const {
	return allReactants;
}

/**
 * This operation returns all reactants in the network with the given name.
 * The list may or may not be ordered and the decision is left to
 * implementers.
 * @param name The reactant or compound reactant name
 * @return The list of all of the reactants in the network or null if the
 * name is invalid.
 */
std::vector<Reactant *> PSIClusterReactionNetwork::getAll(
		std::string name) const {

	// Local Declarations
	std::vector<Reactant *> reactants;

	// Only pull the reactants if the name is valid
	if (name == heType || name == vType || name == iType || name == "HeV"
			|| name == "HeI") {
		std::shared_ptr < std::vector<std::shared_ptr<Reactant>>
				> storedReactants = clusterTypeMap.at(name);
		int vecSize = storedReactants->size();
		for (int i = 0; i < vecSize; i++) {
			reactants.push_back(storedReactants->at(i).get());
		}
	}

	return reactants;
}

/**
 * This operation adds a reactant or a compound reactant to the network.
 * @param reactant The reactant that should be added to the network.
 */
void PSIClusterReactionNetwork::add(std::shared_ptr<Reactant> reactant) {

	// Local Declarations
	int numHe = 0, numV = 0, numI = 0;
	bool isMixed = false;
	std::string numClusterKey, clusterSizeKey, name;

	// Only add a complete reactant
	if (reactant != NULL) {
		// Check the type... but call it name... whatever.
		name = reactant->getName();
		// Get the composition
		auto composition = reactant->getComposition();
		// Get the species sizes
		numHe = composition.at(heType);
		numV = composition.at(vType);
		numI = composition.at(iType);
		// Determine if the cluster is a compound. If there is more than one
		// type, then the check below will sum to greater than one and we know
		// that we have a mixed cluster.
		isMixed = ((numHe > 0) + (numV > 0) + (numI > 0)) > 1;
		// Only add the element if we don't already have it
		// Add the compound or regular reactant.
		if (isMixed && mixedSpeciesMap.count(composition) == 0) {
			// Put the compound in its map
			mixedSpeciesMap[composition] = std::dynamic_pointer_cast
					< PSICluster > (reactant);
			// Figure out whether we have HeV or HeI and set the keys
			if (numV > 0) {
				numClusterKey = "numHeVClusters";
				clusterSizeKey = "maxHeVClusterSize";
			} else {
				numClusterKey = "numHeIClusters";
				clusterSizeKey = "maxHeIClusterSize";
			}
		} else if (!isMixed && singleSpeciesMap.count(composition) == 0) {
			/// Put the reactant in its map
			singleSpeciesMap[composition] = std::dynamic_pointer_cast
					< PSICluster > (reactant);
			// Figure out whether we have He, V or I and set the keys
			if (numHe > 0) {
				numClusterKey = "numHeClusters";
				clusterSizeKey = "maxHeClusterSize";
			} else if (numV > 0) {
				numClusterKey = "numVClusters";
				clusterSizeKey = "maxVClusterSize";
			} else {
				numClusterKey = "numIClusters";
				clusterSizeKey = "maxIClusterSize";
			}
		} else {
			std::stringstream errStream;
			errStream << "PSIClusterReactionNetwork Message: "
					<< "Duplicate Reactant (He=" << numHe << ",V=" << numV
					<< ",I=" << numI << ") not added!" << std::endl;
			throw errStream.str();
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
		// Update the size
		++networkSize;
		// Set the id for this cluster
		reactant->setId(networkSize);
		// Get the vector for this reactant from the type map
		auto clusters = clusterTypeMap[reactant->getType()];
		clusters->push_back(reactant);
		// Add the pointer to the list of all clusters
		allReactants->push_back(reactant.get());
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
void PSIClusterReactionNetwork::setProperty(std::string key,
		std::string value) {

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

/**
 * This operation returns the size or number of reactants in the network.
 * @return The number of reactants in the network
 */int PSIClusterReactionNetwork::size() {
	return networkSize;
}
