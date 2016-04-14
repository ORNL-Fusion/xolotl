#include "PSIClusterReactionNetwork.h"
#include "PSICluster.h"
#include <xolotlPerf.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <Constants.h>

using namespace xolotlCore;

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
		ReactionNetwork(), temperature(0.0) {
	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork(registry), temperature(0.0) {
	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork(
		const PSIClusterReactionNetwork &other) :
		ReactionNetwork(other),
		names(other.getNames()),
		compoundNames(other.getCompoundNames()),
		temperature(other.getTemperature()) {
	// The size and ids do not need to be copied. They will be fixed when the
	// reactants are added.

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
	for (unsigned int i = 0; i < reactants.size(); i++) {
		add(reactants[i]->clone());
	}

	return;
}

void PSIClusterReactionNetwork::setTemperature(double temp) {
	// Set the temperature
	temperature = temp;

	// Update the temperature for all of the clusters
	int networkSize = size();
	for (int i = 0; i < networkSize; i++) {
		// This part will set the temperature in each reactant
		// and recompute the diffusion coefficient
		allReactants->at(i)->setTemperature(temp);
	}
	for (int i = 0; i < networkSize; i++) {
		// Now that the diffusion coefficients of all the reactants
		// are updated, the reaction and dissociation rates can be
		// recomputed
		auto cluster = (PSICluster *) allReactants->at(i);
		cluster->computeRateConstants();
	}

	return;
}

double PSIClusterReactionNetwork::getTemperature() const {
	return temperature;
}

Reactant * PSIClusterReactionNetwork::get(const std::string& type,
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

Reactant * PSIClusterReactionNetwork::getCompound(const std::string& type,
		const std::vector<int>& sizes) const {
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

const std::shared_ptr<std::vector<Reactant *>> & PSIClusterReactionNetwork::getAll() const {
	return allReactants;
}

std::vector<Reactant *> PSIClusterReactionNetwork::getAll(
		const std::string& name) const {
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

void PSIClusterReactionNetwork::add(std::shared_ptr<Reactant> reactant) {
	// Local Declarations
	int numHe = 0, numV = 0, numI = 0;
	bool isMixed = false;
	std::string numClusterKey, clusterSizeKey;

	// Only add a complete reactant
	if (reactant != NULL) {
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
		}
		else if (!isMixed && singleSpeciesMap.count(composition) == 0) {
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
		}
		else {
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

void PSIClusterReactionNetwork::reinitializeConnectivities() {
	// Loop on all the reactants to reset their connectivities
	PSICluster * cluster;
	for (auto it = allReactants->begin(); it != allReactants->end(); ++it) {
		cluster = (PSICluster *) *it;
		cluster->resetConnectivities();
	}

	return;
}

const std::vector<std::string> & PSIClusterReactionNetwork::getNames() const {
	return names;
}

const std::vector<std::string> & PSIClusterReactionNetwork::getCompoundNames() const {
	return compoundNames;
}

const std::map<std::string, std::string> & PSIClusterReactionNetwork::getProperties() {
	return *properties;
}

void PSIClusterReactionNetwork::setProperty(const std::string& key,
		const std::string& value) {
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

int PSIClusterReactionNetwork::size() {
	return networkSize;
}


