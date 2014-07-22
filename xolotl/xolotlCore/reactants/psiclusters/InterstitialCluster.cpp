#include "InterstitialCluster.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

InterstitialCluster::InterstitialCluster(int nI, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(nI, registry) {

	// Update the composition map
	compositionMap["I"] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "I_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "I";

	// Compute the reaction radius
	double EightPi = 8.0 * xolotlCore::pi;
	double aCubed = pow(xolotlCore::latticeConstant, 3.0);
	double termOne = 1.15 * (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant;
	double termTwo = pow((3.0 / EightPi) * aCubed * size, (1.0 / 3.0));
	double termThree = pow((3.0 / EightPi) * aCubed, (1.0 / 3.0));
	reactionRadius = termOne + termTwo - termThree;

}

InterstitialCluster::~InterstitialCluster() {
}

std::shared_ptr<Reactant> InterstitialCluster::clone() {
	std::shared_ptr<Reactant> reactant(new InterstitialCluster(*this));
	return reactant;
}

void InterstitialCluster::createReactionConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species reaction
	PSICluster::createReactionConnectivity();

	// This cluster is always I_a

	// Helium-Interstitial clustering
	// He_b + I_a --> (He_b)(I_a)
	// Get all the He clusters from the network
	auto reactants = network->getAll(heType);
	// combineClusters handles He combining with I to form HeI
	combineClusters(reactants, heIType);

	// Single Interstitial absorption by HeI clusters
	// (He_b)(I_c) + I_a --> (He_b)[I_(c+a)]
	// for a = 1
	// Only if the size of this cluster is 1
	if (size == 1) {
		// Get all the HeI clusters from the network
		reactants = network->getAll(heIType);
		// combineClusters handles HeI combining with I to form HeI
		combineClusters(reactants, heIType);
	}

	// Vacancy-Interstitial annihilation
	// I_a + V_b
	//        --> I_(a-b), if a > b
	//        --> V_(b-a), if a < b
	//        --> 0, if a = b
	// Get all the V clusters from the network
	reactants = network->getAll(vType);
	// fillVWithI handles this reaction
	fillVWithI(vType, reactants);

	// Vacancy-Interstitial annihilation producing this cluster
	// I_(a+b) + V_b --> I_a
	// All the V clusters are already in reactants
	auto reactantsSize = reactants.size();
	for (int i = 0; i < reactantsSize; i++) {
		auto firstReactant = (PSICluster *) reactants[i];
		// Get the interstitial cluster that is bigger than the vacancy
		// and can form this cluster.
		auto secondReactant = (PSICluster *) network->get(typeName,firstReactant->getSize() + size);
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

	// Vacancy reduction by Interstitial absorption in HeV clusters
	// (He_b)(V_c) + (I_a) --> (He_b)[V_(c-a)]
	// Get all the HeV clusters from the network
	reactants = network->getAll(heVType);
	// replaceInCompound handles this reaction
	replaceInCompound(reactants, vType, iType);

	return;
}

void InterstitialCluster::createDissociationConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species dissociation
	PSICluster::createDissociationConnectivity();

	// This cluster is always I_a

	// Specific case for the single species cluster
	if (size == 1) {
		// I dissociation of HeI cluster is handled here
		// (He_c)(I_b) --> (He_c)[I_(b-a)] + I_a
		// for a = 1
		// Get all the HeI clusters of the network
		auto allHeIReactants = network->getAll(heIType);
		for (int i = 0; i < allHeIReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeIReactants[i];

			// (He_c)(I_b) is the dissociating one, (He_c)[I_(b-a)] is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();
			std::vector<int> compositionVec = { comp[heType], comp[vType],
					comp[iType] - 1 };
			auto smallerReactant = network->getCompound(heIType, compositionVec);
			dissociateCluster(cluster, smallerReactant);
		}

		// Trap mutation of HeV cluster is handled here
		// (He_c)(V_b) --> He_(c)[V_(b+a)] + I_a
		// for a = 1
		// Get all the HeV clusters of the network
		auto allHeVReactants = network->getAll(heVType);
		for (int i = 0; i < allHeVReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeVReactants[i];

			// (He_c)(V_b) is the dissociating one, (He_c)[V_(b+a)] is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();
			std::vector<int> compositionVec = { comp[heType], comp[vType] + 1,
					comp[iType] };
			auto biggerReactant = network->getCompound(heVType, compositionVec);
			dissociateCluster(cluster, biggerReactant);
		}
	}

	return;
}
