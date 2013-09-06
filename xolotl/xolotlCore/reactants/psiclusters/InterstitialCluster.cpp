#include "InterstitialCluster.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

InterstitialCluster::InterstitialCluster(int nI) :
		PSICluster(nI) {
	// Set the reactant name appropriately
	name = "I";
}

InterstitialCluster::~InterstitialCluster() {
}

std::shared_ptr<Reactant> InterstitialCluster::clone() {
	std::shared_ptr<Reactant> reactant(new InterstitialCluster(*this));
	return reactant;
}

void InterstitialCluster::createReactionConnectivity() {

	// Local Declarations - Note the reference to the properties map
	auto props = network->getProperties();
	int numHe, indexOther, networkSize = network->size();
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxIClusterSize = std::stoi(props["maxIClusterSize"]);
	int maxHeIClusterSize = std::stoi(props["maxHeIClusterSize"]);
	int numHeVClusters = std::stoi(props["numHeVClusters"]);
	int numHeIClusters = std::stoi(props["numHeIClusters"]);
	int numIClusters = std::stoi(props["numIClusters"]);
	std::map<std::string, int> composition;
	int totalSize = 1, firstSize = 0, secondSize = 0;
	int firstIndex = -1, secondIndex = -1, reactantVecSize = 0;
	std::shared_ptr<Reactant> firstReactant, secondReactant;
	std::shared_ptr<PSICluster> psiCluster;

	/*
	 * This section fills the array of reacting pairs that combine to produce
	 * this cluster. The only reactions that produce I clusters are those I
	 * clusters that are smaller than this.size. Each cluster i combines with
	 * a second cluster of size this.size - i.size.
	 *
	 * Total size starts with a value of one so that clusters of size one are
	 * not considered in this loop.
	 */
	while (totalSize < size) {
		// Increment the base sizes
		++firstSize;
		secondSize = size - firstSize;
		// Get the first and second reactants for the reaction
		// first + second = this.
		firstReactant = network->get("I", firstSize);
		secondReactant = network->get("I", secondSize);
		// Create a ReactingPair with the two reactants
		if (firstReactant && secondReactant) {
			ReactingPair pair;
			pair.first = std::dynamic_pointer_cast<PSICluster>(firstReactant);
			pair.second = std::dynamic_pointer_cast<PSICluster>(secondReactant);
			// Add the pair to the list
			reactingPairs.push_back(pair);
		}
		// Update the total size. Do not delete this or you'll have an infinite
		// loop!
		totalSize = firstSize + secondSize;
	}

	/* ----- A*I + B*I → (A+B)*I -----
	 *	Interstitial absorption
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	auto reactants = network->getAll("I");
	reactantVecSize = reactants->size();
	for (int i = 0; i < reactantVecSize; i++) {
		// Get the reactant, its composition and id
		firstReactant = reactants->at(i);
		composition = firstReactant->getComposition();
		indexOther = network->getReactantId(*firstReactant) - 1;
		// React if the size of the product is valid
		if ((size + composition["I"] <= maxIClusterSize)) {
			reactionConnectivity[indexOther] = 1;
			combiningReactants.push_back(firstReactant);
		}
	}

	/* ----- (A*He) + (B*I) --> (A*He)*(B*I)
	 * Interstitials can interact with clusters of He to form HeI clusters.
	 * They cannot cluster with He clusters that are so large that the
	 * combination of the two would produce an HeI cluster above the
	 * maximum size.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	reactants = network->getAll("He");
	reactantVecSize = reactants->size();
	for (int i = 0; i < reactantVecSize; i++) {
		// Get the reactant, its composition and id
		firstReactant = reactants->at(i);
		composition = firstReactant->getComposition();
		indexOther = network->getReactantId(*firstReactant) - 1;
		// React if the size of the product is valid
		if ((size + composition["He"] <= maxHeIClusterSize)) {
			reactionConnectivity[indexOther] = 1;
			combiningReactants.push_back(firstReactant);
		}
	}

	/* ----- A*I + B*V -----
	 * → (A-B)*I, if A > B
	 * → (B-I)*V, if A < B
	 * → 0, if A = B
	 * Interstitial-Vacancy Annihilation
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	reactants = network->getAll("V");
	reactantVecSize = reactants->size();
	for (int i = 0; i < reactantVecSize; i++) {
		// Get the reactant and its id
		firstReactant = reactants->at(i);
		indexOther = network->getReactantId(*firstReactant) - 1;
		// Always interact with vacancies
		reactionConnectivity[indexOther] = 1;
		combiningReactants.push_back(firstReactant);
	}

	/* ----- (A*He)(B*V) + (C*I) --> (A*He)[(B-C)V] -----
	 * Interstitials interact with all mixed-species clusters by
	 * annihilating vacancies.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	if (numHeIClusters > 0) {
		reactants = network->getAll("HeV");
		reactantVecSize = reactants->size();
		for (int i = 0; i < reactantVecSize; i++) {
			// Get the reactant and its id
			firstReactant = reactants->at(i);
			indexOther = network->getReactantId(*firstReactant) - 1;
			// Always interact with HeV
			reactionConnectivity[indexOther] = 1;
			combiningReactants.push_back(firstReactant);
		}
	}

	/* ----- (A*He)*(B*I) + I --> (A*He)*(B + 1)*I -----
	 * Single interstitial absorption by a HeI cluster under the condition
	 * that (x + y + 1) <= maxSize
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	if (size == 1 && numHeVClusters > 0) {
		reactants = network->getAll("HeI");
		reactantVecSize = reactants->size();
		for (int i = 0; i < reactantVecSize; i++) {
			// Get the reactant, and its id
			firstReactant = reactants->at(i);
			indexOther = network->getReactantId(*firstReactant) - 1;
			// React if the size of the product is valid
			psiCluster = std::dynamic_pointer_cast<PSICluster>(firstReactant);
			if ((size + psiCluster->getSize() <= maxHeIClusterSize)) {
				reactionConnectivity[indexOther] = 1;
				combiningReactants.push_back(firstReactant);
			}
		}
	}

	return;
}

void InterstitialCluster::createDissociationConnectivity() {

	// Local Declarations
	int nReactants = network->size(), id = 0;
	std::map<std::string, int> clusterMap;
	std::shared_ptr<Reactant> reactant;
	auto props = network->getProperties();

	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	int maxHeVClusterSize = std::stoi(props["maxHeVClusterSize"]);
	int numHeVClusters = std::stoi(props["numHeVClusters"]);
	int numHeIClusters = std::stoi(props["numHeIClusters"]);
	int numIClusters = std::stoi(props["numIClusters"]);

	// Interstitial dissociation, get a vacancy with size = size - 1
	reactant = network->get("I", size - 1);
	if (reactant) {
		id = network->getReactantId(*reactant);
		dissociationConnectivity[id] = 1;
		// Single V
		reactant = network->get("I", 1);
		id = network->getReactantId(*reactant);
		dissociationConnectivity[id] = 1;
	}

	return;
}

bool InterstitialCluster::isProductReactant(const Reactant & reactantI,
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

	// We should have this->size interstitials, a
	// total of 0 Helium, and a total of
	// 0 Vacancies
	return ((rI_I + rJ_I) == size) && ((rI_He + rJ_He) == 0)
			&& ((rI_V + rJ_V) == 0);
}

std::map<std::string, int> InterstitialCluster::getClusterMap() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = 0;
	clusterMap["V"] = 0;
	clusterMap["I"] = size;

	// Return it
	return clusterMap;
}

std::map<std::string, int> InterstitialCluster::getComposition() const {

	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = 0;
	clusterMap["V"] = 0;
	clusterMap["I"] = size;

	// Return it
	return clusterMap;
}

double InterstitialCluster::getReactionRadius() {

	double EightPi = 8.0 * xolotlCore::pi;
	double aCubed = pow(xolotlCore::latticeConstant, 3.0);
	double termOne = 1.15 * (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant;
	double termTwo = pow((3.0 / EightPi) * aCubed * size, (1.0 / 3.0));
	double termThree = pow((3.0 / EightPi) * aCubed, (1.0 / 3.0));

	return termOne + termTwo - termThree;
}
