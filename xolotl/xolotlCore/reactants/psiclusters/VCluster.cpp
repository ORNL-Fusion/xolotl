// Includes
#include "VCluster.h"
#include <iostream>
#include <Constants.h>
#include <PSIClusterReactionNetwork.h>

using namespace xolotlCore;

VCluster::VCluster(int nV,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(registry) {
	// Set the size
	size = nV;
	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "V_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = vType;

	// Update the composition map
	compositionMap[vType] = size;

	// Compute the reaction radius
	// It is the same formula for HeV clusters
	reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::tungstenLatticeConstant
			+ pow(
					(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0) * size)
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
			- pow(
					(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0))
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0));
}

void VCluster::createReactionConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species reaction
	PSICluster::createReactionConnectivity();

	// This cluster is always V_a

	// Helium-Vacancy clustering
	// He_b + V_a --> (He_b)(V_a)
	// Get all the He clusters from the network
	auto reactants = network->getAll(heType);
	// combineClusters handles He combining with V to form HeV
	combineClusters(reactants, heVType);

	// Single Vacancy absorption by HeV clusters
	// (He_c)(V_b) + V_a --> (He_c)[V_(b+a)]
	// for a = 1
	// Only if the size of this cluster is 1
	if (size == 1) {
		// Get all the HeV clusters from the network
		reactants = network->getAll(heVType);
		// combineClusters handles HeV combining with V to form HeV
		combineClusters(reactants, heVType);
	}

	// Vacancy-Interstitial annihilation
	// I_a + V_b
	//        --> I_(a-b), if a > b
	//        --> V_(b-a), if a < b
	//        --> 0, if a = b
	// Get all the I clusters from the network
	reactants = network->getAll(iType);
	// fillVWithI handles this reaction
	fillVWithI(reactants);

	// Vacancy-Interstitial annihilation producing this cluster
	// I_b + V_(a+b) --> V_a
	// All the I clusters are already in reactants
	int reactantsSize = reactants.size();
	for (int i = 0; i < reactantsSize; i++) {
		auto firstReactant = (PSICluster *) reactants[i];
		// Get the vacancy cluster that is bigger than the interstitial
		// and can form this cluster.
		auto secondReactant = (PSICluster *) network->get(typeName,
				firstReactant->getSize() + size);
		// Add to the reacting pairs if the second reactant exists
		if (secondReactant) {
			// Create the pair
			// The reaction constant will be computed later, it is set to 0.0 for now
			ClusterPair pair(firstReactant, secondReactant, 0.0);
			// Add the pair to the list
			reactingPairs.push_back(pair);
			// Setup the connectivity array
			int Id = firstReactant->getId();
			setReactionConnectivity(Id);
			Id = secondReactant->getId();
			setReactionConnectivity(Id);
		}
	}

	// Interstitial reduction by Vacancy absorption in HeI clusters
	// (He_c)*(I_b) + (V_a) --> (He_c)*[I_(b-a)]
	// Get all the HeI clusters from the network
	reactants = network->getAll(heIType);
	// replaceInCompound handles this reaction
	replaceInCompound(reactants, iType);

	return;
}

void VCluster::createDissociationConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species dissociation
	PSICluster::createDissociationConnectivity();

	// This cluster is always V_a

	// He Dissociation
	// (He_1)(V_a) --> V_a + He
	// Get the cluster with one more helium
	std::vector<int> compositionVec = { 1, size, 0 };
	auto heVClusterMoreHe = (PSICluster *) network->getCompound(heVType,
			compositionVec);
	// Get the single helium cluster
	auto singleCluster = (PSICluster *) network->get(heType, 1);
	// Here it is important that heVClusterMoreHe is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heVClusterMoreHe, singleCluster);

	// Specific case for the single species cluster
	if (size == 1) {
		// V dissociation of HeV cluster is handled here
		// (He_c)(V_b) --> (He_c)[V_(b-a)] + V_a
		// for a = 1
		// Get all the HeV clusters of the network
		auto allHeVReactants = network->getAll(heVType);
		for (unsigned int i = 0; i < allHeVReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeVReactants[i];

			// (He_c)(V_b) is the dissociating one, (He_c)[V_(b-a)] is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();

			// Skip He_1V_1 because it was counted in the He dissociation
			if (comp[heType] == 1 && comp[vType] == 1)
				continue;

			std::vector<int> compositionVec = { comp[heType], comp[vType]
					- size, 0 };
			auto smallerReactant = (PSICluster *) network->getCompound(heVType,
					compositionVec);
			// Special case for numV = 1
			if (comp[vType] == 1) {
				smallerReactant = (PSICluster *) network->get(heType,
						comp[heType]);
			}
			dissociateCluster(cluster, smallerReactant);
		}
	}

	return;
}
