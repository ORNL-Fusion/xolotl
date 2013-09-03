// Includes
#include "HeVCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>

using namespace xolotlCore;

HeVCluster::HeVCluster(int numHe, int numV) :
		PSICluster(1), numHe(numHe), numV(numV) {

	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = numHe + numV;

	// Set the reactant name appropriately
	name = "HeV";
}

HeVCluster::HeVCluster(const HeVCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numV = other.numV;
}

HeVCluster::~HeVCluster() {
}

std::shared_ptr<Reactant> HeVCluster::clone() {
	std::shared_ptr<Reactant> reactant(new HeVCluster(*this));
	return reactant;
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

	// Local Declarations
	auto psiNetwork = std::dynamic_pointer_cast < PSIClusterReactionNetwork
			> (network);
	auto props = psiNetwork->getProperties();
	int networkSize = psiNetwork->size(), index = 0;
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	int maxIClusterSize = std::stoi(props["maxIClusterSize"]);
	int maxHeVClusterSize = std::stoi(props["maxHeVClusterSize"]);
	int maxHeIClusterSize = std::stoi(props["maxHeIClusterSize"]);
	int numHeClusters = std::stoi(props["numHeClusters"]);
	int numVClusters = std::stoi(props["numVClusters"]);
	int numIClusters = std::stoi(props["numIClusters"]);
	int numSingleSpeciesClusters = numHeClusters + numVClusters + numIClusters;
	std::shared_ptr<Reactant> firstReactant, secondReactant;
	std::shared_ptr<PSICluster> heCluster, vCluster, iCluster;
	std::vector<int> firstComposition, secondComposition, speciesMap;

	std::cout << numHe << " " << numV << " " << std::endl;

	/* ----- (A*He)(B*V) + (C*He) --> [(A+C)He]*(B*V) -----
	 * Helium absorption by HeV clusters that results
	 * in the production of this cluster.
	 */
	for (int z = 1; z <= maxHeClusterSize; z++) {
		// Get the first reactant
		firstComposition = psiNetwork->getCompositionVector(numHe - z, numV, 0);
		firstReactant = psiNetwork->getCompound("HeV", firstComposition);
		// Get the second reactant
		secondReactant = psiNetwork->get("He", z);
		// Create the Reacting Pair
		ReactingPair pair;
		pair.first = std::dynamic_pointer_cast < PSICluster > (firstReactant);
		pair.second = std::dynamic_pointer_cast < PSICluster > (secondReactant);
		// Add the pair to the list
		reactingPairs.push_back(pair);
	}

	/* ---- (AHe)*(BV) + (CHe) --> [(A + C)He]*(BV) ----
	 * HeV clusters can absorb helium clusters so long as they do not cross
	 * the max size limit.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	auto reactants = psiNetwork->getAll("He");
	int numReactants = reactants->size();
	for (int i = 0; i < numReactants; i++) {
		heCluster = std::dynamic_pointer_cast < PSICluster > (reactants->at(i));
		// React with it if the sizes are compatible
		if (heCluster->getSize() + size <= maxHeVClusterSize) {
			index = psiNetwork->getReactantId(*heCluster) - 1;
			reactionConnectivity[i] = 1;
			combiningReactants.push_back(heCluster);
		}
	}

	/* ----- (A*He)(B*V) + V --> (A*He)[(B+1)*V] -----
	 * Vacancies are also produced by single-vacancy absorption by an HeV
	 * cluster. In this case, (A*He)[(B-1)*V] produces the current cluster.
	 */
	firstComposition = psiNetwork->getCompositionVector(numHe, numV - 1, 0);
	firstReactant = psiNetwork->getCompound("HeV", firstComposition);
	secondReactant = psiNetwork->get("V", 1);
	// Create the Reacting Pair
	ReactingPair pair;
	pair.first = std::dynamic_pointer_cast < PSICluster > (firstReactant);
	pair.second = std::dynamic_pointer_cast < PSICluster > (secondReactant);
	// Add the pair to the list
	reactingPairs.push_back(pair);

	/* ----- (A*He)(B*V) + C*I --> (A*He)[(B-C)V] -----
	 * Interstitial absorption by an HeV cluster produces an HeV cluster of
	 * size B-C vacancies smaller.
	 */
	for (int z = 1; z <= maxIClusterSize; z++) {
		// Get the first reactant
		firstComposition = psiNetwork->getCompositionVector(numHe, numV + z, 0);
		firstReactant = psiNetwork->getCompound("HeV", firstComposition);
		// Get the second reactant
		secondReactant = psiNetwork->get("I", z);
		// Create the Reacting Pair
		ReactingPair pair;
		pair.first = std::dynamic_pointer_cast < PSICluster > (firstReactant);
		pair.second = std::dynamic_pointer_cast < PSICluster > (secondReactant);
		// Add the pair to the list
		reactingPairs.push_back(pair);
	}

	/* ----- (AHe)*(BV) + V --> (AHe)*(B + 1)V -----
	 * HeV clusters can absorb single vacancies.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	if (numHe + numV + 1 <= maxHeVClusterSize) {
		// Get the HeV cluster that is one He bigger than us.
		firstComposition = psiNetwork->getCompositionVector(numHe + 1, numV, 0);
		firstReactant = psiNetwork->getCompound("HeV", firstComposition);
		// Add it to the list if it exists
		if (firstReactant) {
			index = psiNetwork->getReactantId(*firstReactant) - 1;
			reactionConnectivity.at(index) = 1;
			combiningReactants.push_back(firstReactant);
		}
		// Get the HeV cluster that is on V bigger than us.
		secondComposition = psiNetwork->getCompositionVector(numHe, numV + 1,
				0);
		secondReactant = psiNetwork->getCompound("HeV", secondComposition);
		// Add it to the list if it exists
		if (secondReactant) {
			index = psiNetwork->getReactantId(*secondReactant) - 1;
			reactionConnectivity.at(index) = 1;
			combiningReactants.push_back(secondReactant);
		}
	}

	/* ----- (AHe)*(BV) + CI  --> (AHe)*(B - C)V -----
	 * Helium-vacancy clusters lose vacancies when they interact with
	 * interstitial clusters.
	 *
	 * We assume that the HeV and interstitial cluster can only
	 * interact if they would produce another HeV cluster, not single He.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	reactants = psiNetwork->getAll("I");
	numReactants = reactants->size();
	for (int i = 0; i < numReactants; i++) {
		iCluster = std::dynamic_pointer_cast < PSICluster > (reactants->at(i));
		// Only add it if it exists!
		if (iCluster->getSize() - numV >= 1) {
			index = psiNetwork->getReactantId(*iCluster) - 1;
			reactionConnectivity.at(index) = 1;
			combiningReactants.push_back(firstReactant);
		}
	}

	return;
}

void HeVCluster::createDissociationConnectivity() {

	// Local Declarations
	auto psiNetwork = std::dynamic_pointer_cast < PSIClusterReactionNetwork
			> (network);
	auto props = psiNetwork->getProperties();
	int index = 0;
	std::vector<int> composition;
	std::shared_ptr<Reactant> singleCluster, otherMixedCluster;

	// He Dissociation, get the [(numHe-1)*He]V and He
	composition = psiNetwork->getCompositionVector(numHe - 1, numV, 0);
	otherMixedCluster = psiNetwork->getCompound("HeV", composition);
	singleCluster = psiNetwork->get("He", 1);
	if (singleCluster && otherMixedCluster) {
		// Add the HeV
		index = psiNetwork->getReactantId(*otherMixedCluster) - 1;
		dissociationConnectivity[index] = 1;
		// Add the He
		index = psiNetwork->getReactantId(*singleCluster) - 1;
		dissociationConnectivity[index] = 1;
	}

	// Vacancy Dissociation, get He[(numV-1)*V] and V
	composition = psiNetwork->getCompositionVector(numHe, numV - 1, 0);
	otherMixedCluster = psiNetwork->getCompound("HeV", composition);
	singleCluster = psiNetwork->get("V", 1);
	if (singleCluster && otherMixedCluster) {
		// Handle the HeV
		int index = psiNetwork->getReactantId(*otherMixedCluster) - 1;
		dissociationConnectivity[index] = 1;
		// Handle the V
		index = psiNetwork->getReactantId(*singleCluster) - 1;
		dissociationConnectivity[index] = 1;
	}

	// Trap mutation, get He[(numV+1)*V] and I
	composition = psiNetwork->getCompositionVector(numHe, numV + 1, 0);
	otherMixedCluster = psiNetwork->getCompound("HeV", composition);
	singleCluster = psiNetwork->get("I", 1);
	if (singleCluster && otherMixedCluster) {
		// Handle the HeV
		int index = psiNetwork->getReactantId(*otherMixedCluster) - 1;
		dissociationConnectivity[index] = 1;
		// Handle the I
		index = psiNetwork->getReactantId(*singleCluster) - 1;
		dissociationConnectivity[index] = 1;
	}

	return;
}

double HeVCluster::getDissociationFlux(double temperature) {

	// Local Declarations
	std::map<std::string, int> composition;
	std::shared_ptr<PSICluster> currentReactant, secondReactant;
	double f4 = 0.0, f3 = 0.0;

	// Calculate the much easier f4 term...
	auto heCluster = std::dynamic_pointer_cast < PSICluster
			> (network->get("He", 1));
	auto vCluster = std::dynamic_pointer_cast < PSICluster
			> (network->get("V", 1));
	auto iCluster = std::dynamic_pointer_cast < PSICluster
			> (network->get("I", 1));
	f4 = calculateDissociationConstant(*this, *iCluster, temperature)
			+ calculateDissociationConstant(*this, *vCluster, temperature)
			+ calculateDissociationConstant(*this, *heCluster, temperature);

	// Loop over all the elements of the dissociation
	// connectivity to find where this mixed species dissociates
	auto reactants = network->getAll();
	int numReactants = dissociationConnectivity.size();
	for (int i = 0; i < numReactants; i++) {
		if (dissociationConnectivity[i] == 1) {
			// Set the current reactant
			currentReactant = std::dynamic_pointer_cast < PSICluster
					> (reactants->at(i));
			// Get the cluster map of this connection
			composition = currentReactant->getComposition();
			// We need to find if this is a Helium dissociation
			if (numHe - composition["He"] == 1 && numV == composition["V"]
					&& composition["I"] == 0) {
				secondReactant = heCluster;
			} else if (numHe == composition["He"] && numV - composition["V"] == 1
					&& composition["V"] == 0) {
				// vacancy dissociation
				secondReactant = vCluster;
			} else if (numHe == composition["He"] && composition["I"] - numV == 1
					&& composition["V"] == 0) {
				// or a trap mutation.
				secondReactant = iCluster;
			}
			// Update the flux calculation
			f3 += calculateDissociationConstant(*currentReactant, *secondReactant,
					temperature) * currentReactant->getConcentration();
		}
	}

	return f3 - f4 * getConcentration();
}

bool HeVCluster::isProductReactant(const Reactant & reactantI,
		const Reactant & reactantJ) {
	// Local Declarations, integers for species number for I, J reactants
	int rI_I = 0, rJ_I = 0, rI_He = 0, rJ_He = 0, rI_V = 0, rJ_V = 0;

	// Get the compositions of the reactants
	auto reactantIMap = reactantI.getComposition();
	auto reactantJMap = reactantJ.getComposition();

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

std::map<std::string, int> HeVCluster::getComposition() const {
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
