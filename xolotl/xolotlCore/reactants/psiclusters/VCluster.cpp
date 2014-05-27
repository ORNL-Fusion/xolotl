// Includes
#include "VCluster.h"
#include <iostream>
#include <Constants.h>
#include <PSIClusterReactionNetwork.h>

using namespace xolotlCore;

VCluster::VCluster(int nV, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(nV, registry) {

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "V_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "V";

	// Update the composition map
	compositionMap["V"] = size;

	// Compute the reaction radius
	// FIXME Not right...
	reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant;

}

VCluster::~VCluster() {
}

std::shared_ptr<Reactant> VCluster::clone() {
	std::shared_ptr<Reactant> reactant(new VCluster(*this));
	return reactant;
}

void VCluster::createReactionConnectivity() {

	// Local Declarations
	auto props = network->getProperties();
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	int maxHeVClusterSize = std::stoi(props["maxHeVClusterSize"]);
	int numHeVClusters = std::stoi(props["numHeVClusters"]);
	int numHeIClusters = std::stoi(props["numHeIClusters"]);
	int firstSize = 0, secondSize = 0;

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(getId());

	/*
	 * This section fills the array of reacting pairs that combine to produce
	 * this cluster. The only reactions that produce V clusters are those V
	 * clusters that are smaller than this.size. Each cluster i combines with
	 * a second cluster of size this.size - i.size.
	 *
	 * Total size starts with a value of one so that clusters of size one are
	 * not considered in this loop.
	 */
	for (firstSize = 1; firstSize <= (int) size/2; firstSize++) {
		secondSize = size - firstSize;
		// Get the first and second reactants for the reaction
		// first + second = this.
		auto firstReactant = (PSICluster *) network->get("V", firstSize);
		auto secondReactant = (PSICluster *) network->get("V", secondSize);
		// Create a ReactingPair with the two reactants
		if (firstReactant && secondReactant) {
			ReactingPair pair(firstReactant,secondReactant);
			// Add the pair to the list
			reactingPairs.push_back(pair);
		}
	}

	/* -----  He_a + V_b --> (He_a)(V_b) -----
	 * Vacancy clusters can interact with any helium cluster so long as the sum
	 * of the number of helium atoms and vacancies does not produce a cluster
	 * with a size greater than the maximum mixed-species cluster size.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	auto reactants = network->getAll("He");
	combineClusters(reactants, maxHeVClusterSize, "HeV");

	/* ----- V_a + V_b --> V_(a+b) -----
	 * This cluster should interact with all other clusters of the same type up
	 * to the max size minus the size of this one to produce larger clusters.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	reactants = network->getAll("V");
	combineClusters(reactants, maxVClusterSize, "V");

	/* ----- I_a + V_b -----
	 * --> I_(a-b), if a > b
	 * --> V_(b-a), if a < b
	 * --> 0, if a = b -----
	 *
	 * Vacancies are always filled by interstitials.
	 */
	reactants = network->getAll("I");
	fillVWithI("I", reactants);
	// Mark the reaction connectivity for the cases where this cluster is
	// produced by the above reaction. This has to be checked for every
	// vacancy.
	auto reactantsSize = reactants.size();
	for (int i = 0; i < reactantsSize; i++) {
		auto firstReactant = (PSICluster *) reactants[i];
		// Get the interstitial cluster that is bigger than the vacancy
		// and can form this cluster. V only results when it is bigger than I.
		auto secondReactant = (PSICluster *) network->get("V",firstReactant->getSize() + size);
		// Update the connectivity
		if (secondReactant) {
			setReactionConnectivity(firstReactant->getId());
			setReactionConnectivity(secondReactant->getId());
		}
	}

	/* ----- (He_a)(V_b) + V --> (He_a)[V_(b+1)] -----
	 * Single vacancies can interact with a mixed-species cluster so long as
	 * the sum of the number of vacancy atoms and the size of the mixed-species
	 * cluster does not exceed the maximum mixed-species cluster size.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	if (size == 1 && numHeVClusters > 0) {
		reactants = network->getAll("HeV");
		combineClusters(reactants, maxHeVClusterSize, "HeV");
	}

	/* ----- (He_a)*(I_b) + (V_c) --> (He_a)*[I_(b-c)] -----
	 * Vacancy absorption by HeI under the condition that y - z >= 1
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	if (numHeIClusters > 0) {
		reactants = network->getAll("HeI");
		replaceInCompound(reactants, "I", "V");
	}

	return;
}
