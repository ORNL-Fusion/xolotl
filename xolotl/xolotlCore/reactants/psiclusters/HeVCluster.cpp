// Includes
#include "HeVCluster.h"
#include <iostream>
#include <Constants.h>

using namespace xolotlCore;

HeVCluster::HeVCluster(int numHe, int numV) :
		PSICluster(1), numHe(numHe), numV(numV) {

	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = numHe + numV;

	// Set the reactant name appropriately
	name = "HeV Cluster";
}

HeVCluster::~HeVCluster() {
}

double HeVCluster::getGenByEm() {
	return 0;
}

double HeVCluster::getAnnByEm() {
	return 0;
}

int HeVCluster::getSpeciesSize(const std::string speciesName) {
	if (speciesName == "He") {
		return numHe;
	} else if (speciesName == "V") {
		return numV;
	} else {
		return 0;
	}
}

void HeVCluster::createReactionConnectivity() {

	// Extract some of the properties from the network
	std::shared_ptr<std::map<std::string, std::string>> properties =
			network->properties;
	int maxMixedClusterSize = std::stoi(properties->at("maxMixedClusterSize"));
	int maxHeClusterSize = std::stoi(properties->at("maxHeClusterSize"));
	int maxIClusterSize = std::stoi(properties->at("maxIClusterSize"));
	int numHeClusters = std::stoi(properties->at("numHeClusters"));
	std::shared_ptr<Reactant> firstReactant, secondReactant;
	std::map<std::string, int> firstReactantMap, secondReactantMap;
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>>reactants =
			network->reactants;
	std::map<std::string, int> speciesMap;
	std::shared_ptr<PSICluster> heCluster;

	std::cout << numHe << " " << numV << " " << std::endl;

	/* ----- (A*He)(B*V) + (C*He) --> [(A+C)He]*(B*V) -----
	 * Fill reacting pairs for helium absorption by mixed clusters that results
	 * in the production of this cluster.
	 */
	for (int z = 1; z <= maxHeClusterSize; z++) {
		// Set the first reactant's map data
		firstReactantMap["He"] = numHe - z;
		firstReactantMap["V"] = numV;
		firstReactantMap["I"] = 0;

		// Set the second's data
		secondReactantMap["He"] = z;
		secondReactantMap["V"] = 0;
		secondReactantMap["I"] = 0;

		int firstIndex = network->toClusterIndex(firstReactantMap);
		int secondIndex = network->toClusterIndex(secondReactantMap);

		if (firstIndex < reactants->size() && secondIndex < reactants->size()) {
			// Get those Reactants from the network
			firstReactant = reactants->at(firstIndex);
			secondReactant = reactants->at(secondIndex);

			// Create the Reacting Pair
			ReactingPair pair;
			pair.first = std::dynamic_pointer_cast < PSICluster
					> (reactants->at(firstIndex));
			pair.second = std::dynamic_pointer_cast < PSICluster
					> (reactants->at(secondIndex));

			// Add the pair to the list
			reactingPairs.push_back(pair);
		}
	}

	// xHe yV + V --> xHe (y+1) V
	// Set the first reactant's map data
	firstReactantMap["He"] = numHe;
	firstReactantMap["V"] = numV - 1;
	firstReactantMap["I"] = 0;

	// Set the second's data
	secondReactantMap["He"] = 0;
	secondReactantMap["V"] = 1;
	secondReactantMap["I"] = 0;

	int firstIndex = network->toClusterIndex(firstReactantMap);
	int secondIndex = network->toClusterIndex(secondReactantMap);

	if (firstIndex < reactants->size() && secondIndex < reactants->size()) {
		// Get those Reactants from the network
		firstReactant = reactants->at(firstIndex);
		secondReactant = reactants->at(secondIndex);

		// Create the Reacting Pair
		ReactingPair pair;
		pair.first = std::dynamic_pointer_cast < PSICluster
				> (reactants->at(firstIndex));
		pair.second = std::dynamic_pointer_cast < PSICluster
				> (reactants->at(secondIndex));

		// Add the pair to the list
		reactingPairs.push_back(pair);
	}

	// (A*He)(B*V) + C*I --> (A*He)[(B-C)V]
	for (int z = 1; z <= maxIClusterSize; z++) {
		// Set the first reactant's map data
		firstReactantMap["He"] = numHe;
		firstReactantMap["V"] = numV + z;
		firstReactantMap["I"] = 0;

		// Set the second's data
		secondReactantMap["He"] = 0;
		secondReactantMap["V"] = 0;
		secondReactantMap["I"] = z;

		// Get those Reactants from the network
		int firstReactantIndex = network->toClusterIndex(firstReactantMap);
		int secondReactantIndex = network->toClusterIndex(secondReactantMap);

		if (firstReactantIndex < reactants->size()
				&& secondReactantIndex < reactants->size()) {
			firstReactant = reactants->at(firstReactantIndex);
			secondReactant = reactants->at(secondReactantIndex);

			// Create the Reacting Pair
			ReactingPair pair;
			pair.first = std::dynamic_pointer_cast < PSICluster
					> (firstReactant);
			pair.second = std::dynamic_pointer_cast < PSICluster
					> (secondReactant);

			// Add the pair to the list
			reactingPairs.push_back(pair);
		}
	}

	/* ---- (AHe)*(BV) + (CHe) --> [(A + C)He]*(BV) ----
	 * HeV clusters can absorb helium clusters so long as they do not cross
	 * the max size.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 *
	 * Start at single species He and count up.
	 */
	speciesMap.clear();
	speciesMap["He"] = 1;
	int heIndex = network->toClusterIndex(speciesMap);
	// Loop over the He clusters
	for (int i = heIndex; i < heIndex + numHeClusters; i++) {
		heCluster = std::dynamic_pointer_cast < PSICluster > (reactants->at(i));
		// React with it if the sizes are compatible.
		if (heCluster->getSize() + size <= maxMixedClusterSize) {
			reactionConnectivity[i] = 1;
			combiningReactants.push_back(reactants->at(i));
		}
	}

	/* ----- (AHe)*(BV) + V --> (AHe)*(B + 1)V -----
	 * HeV clusters can absorb single vacancies.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	if (numHe + numV + 1 <= maxMixedClusterSize) {
		// Get the HeV cluster that is one bigger than us.
		speciesMap.clear();
		speciesMap["He"] = numHe;
		speciesMap["V"] = numV;
		int i = network->toClusterIndex(speciesMap);
		// Set the connectivity and add the reactant to the list.
		reactionConnectivity.at(i) = 1;
		combiningReactants.push_back(reactants->at(i));
	}

	/* ----- (AHe)*(BV) + (CI)  --> (AHe)*(B - C)V -----
	 * Helium-vacancy clusters lose vacancies when they interact with
	 * interstitial clusters.
	 *
	 * We assume that the HeV and interstitial cluster can only
	 * interact if they would produce another HeV cluster, not single He.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	for (int numIOther = 1; numV - numIOther >= 1; numIOther++) {
		// Get the  index
		speciesMap.clear();
		speciesMap["I"] = numIOther;
		int i = network->toClusterIndex(speciesMap);
		reactionConnectivity.at(i) = 1;
		combiningReactants.push_back(reactants->at(i));
	}

	return;
}

void HeVCluster::createDissociationConnectivity() {
// Local Declarations
	std::map<std::string, int> clusterMap;

// Vacancy Dissociation
	clusterMap["He"] = numHe - 1;
	clusterMap["V"] = numV;
	clusterMap["I"] = 0;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["V"] = 0;
	clusterMap["He"] = 1;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;

// Trap Mutation
	clusterMap["V"] = numV + 1;
	clusterMap["He"] = numHe;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["I"] = 1;
	clusterMap["V"] = 0;
	clusterMap["He"] = 0;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;

// Vacancy Dissociation
	clusterMap["He"] = numHe;
	clusterMap["V"] = numV - 1;
	clusterMap["I"] = 0;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["He"] = 0;
	clusterMap["V"] = 1;
	clusterMap["I"] = 0;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
}

double HeVCluster::getDissociationFlux(double temperature) {
// Local Declarations
	std::map<std::string, int> oneHe, oneV, oneI, dissMap;
	std::shared_ptr < std::vector<std::shared_ptr<xolotlCore::Reactant>>
			> reactants;
	std::shared_ptr<Reactant> currentReactant, secondReactant;
	double f4 = 0.0, f3 = 0.0;

// Set the cluster map data for 1 of each species
	oneHe["He"] = 1;
	oneHe["V"] = 0;
	oneHe["I"] = 0;
	oneV["He"] = 0;
	oneV["V"] = 1;
	oneV["I"] = 0;
	oneI["He"] = 0;
	oneI["V"] = 0;
	oneI["I"] = 1;

// Get this PSICluster or subclasses' cluster map
	std::map<std::string, int> thisMap = getClusterMap();

// Get the various indices
	int thisIndex = network->toClusterIndex(thisMap);
	int oneIIndex = network->toClusterIndex(oneI);
	int oneVIndex = network->toClusterIndex(oneV);
	int oneHeIndex = network->toClusterIndex(oneHe);

// Calculate the much easier f4 term...
	reactants = network->reactants;
	f4 = calculateDissociationConstant(reactants->at(thisIndex),
			reactants->at(oneIIndex), temperature)
			+ calculateDissociationConstant(reactants->at(thisIndex),
					reactants->at(oneVIndex), temperature)
			+ calculateDissociationConstant(reactants->at(thisIndex),
					reactants->at(oneHeIndex), temperature);

// Loop over all the elements of the dissociation
// connectivity to find where this mixed species dissociates...
	for (int i = 0; i < dissociationConnectivity.size(); i++) {
		if (dissociationConnectivity[i] == 1) {
			// Set the current reactant
			currentReactant = reactants->at(i);
			// Get the cluster map of this connection
			dissMap = network->toClusterMap(i);
			// We need to find if this is a Helium dissociation,
			// Vacancy dissociation, or a trap mutation.
			if (numHe - dissMap["He"] == 1 && numV == dissMap["I"]
					&& dissMap["V"] == 0) {
				secondReactant = reactants->at(oneHeIndex);
			} else if (numHe == dissMap["He"] && numV - dissMap["V"] == 1
					&& dissMap["V"] == 0) {
				secondReactant = reactants->at(oneVIndex);
			} else if (numHe == dissMap["He"] && dissMap["I"] - numV == 1
					&& dissMap["V"] == 0) {
				secondReactant = reactants->at(oneIIndex);
			}
			// Update the flux calculation
			f3 += calculateDissociationConstant(currentReactant, secondReactant,
					temperature) * currentReactant->getConcentration();
		}
	}

	return f3 - f4 * getConcentration();
}

bool HeVCluster::isProductReactant(int reactantI, int reactantJ) {
// Local Declarations, integers for species number for I, J reactants
	int rI_I = 0, rJ_I = 0, rI_He = 0, rJ_He = 0, rI_V = 0, rJ_V = 0;

// Get the ClusterMap corresponding to
// the given reactants
	std::map<std::string, int> reactantIMap = network->toClusterMap(reactantI);
	std::map<std::string, int> reactantJMap = network->toClusterMap(reactantJ);

// Grab the numbers for each species
// from each Reactant
	rI_I = reactantIMap["I"];
	rJ_I = reactantJMap["I"];
	rI_He = reactantIMap["He"];
	rJ_He = reactantJMap["He"];
	rI_V = reactantIMap["V"];
	rJ_V = reactantJMap["V"];

// We should have no interstitials, a
// total of numHe Helium, and a total of
// numV Vacancies
	return ((rI_I + rJ_I) == 0) && ((rI_He + rJ_He) == numHe)
			&& ((rI_V + rJ_V) == numV);
}

std::map<std::string, int> HeVCluster::getClusterMap() {
// Local Declarations
	std::map<std::string, int> clusterMap;

// Set the number of each species
	clusterMap["He"] = numHe;
	clusterMap["V"] = numV;
	clusterMap["I"] = 0;

// Return it
	return clusterMap;
}

double HeVCluster::getReactionRadius() {
	return (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant
			+ pow(
					(3.0 * pow(xolotlCore::latticeConstant, 3.0) * numV)
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
			- pow(
					(3.0 * pow(xolotlCore::latticeConstant, 3.0))
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0));
}
