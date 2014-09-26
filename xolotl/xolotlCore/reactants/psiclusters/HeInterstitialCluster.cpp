// Includes
#include "HeInterstitialCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>

using namespace xolotlCore;

HeInterstitialCluster::HeInterstitialCluster(int numHelium, int numInterstitial,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(1, registry), numHe(numHelium), numI(numInterstitial) {

	// Set the cluster size as the sum of
	// the number of Helium and Interstitials
	size = numHe + numI;

	// Update the composition map
	compositionMap[heType] = numHe;
	compositionMap[iType] = numI;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << numHe << "I_" << numI;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "HeI";

	// Compute the reaction radius
	reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant
			+ pow(
					(3.0 * pow(xolotlCore::latticeConstant, 3.0) * numI)
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
			- pow(
					(3.0 * pow(xolotlCore::latticeConstant, 3.0))
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0));

	return;
}

HeInterstitialCluster::HeInterstitialCluster(const HeInterstitialCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numI = other.numI;

	return;
}

HeInterstitialCluster::~HeInterstitialCluster() {
}

std::shared_ptr<Reactant> HeInterstitialCluster::clone() {
	std::shared_ptr<Reactant> reactant = std::make_shared<HeInterstitialCluster>(*this);

	return reactant;
}

void HeInterstitialCluster::replaceInCompound(std::vector<Reactant *> & reactants,
		std::string oldComponentName, std::string newComponentName) {
	// Local Declarations
	std::map<std::string, int> myComp = getComposition(),
			productReactantComp;
	int myComponentNumber = myComp[oldComponentName];
	int numReactants = reactants.size();
	int secondId = 0, productId = 0;

	// Loop over all of the extra reactants in this reaction and handle the replacement
	for (int i = 0; i < numReactants; i++) {
		// Get the second reactant and its size
		auto secondReactant = (PSICluster *) reactants[i];
		auto secondReactantSize = secondReactant->getSize();
		// Create the composition vector
		productReactantComp = myComp;
		// Updated the modified components
		productReactantComp[oldComponentName] =
				myComponentNumber - secondReactantSize;
		// Create the composition vector -- FIXME! This should be general!
		std::vector<int> productCompositionVector = { productReactantComp[heType],
				productReactantComp[vType], productReactantComp[iType] };
		// Get the product of the same type as the second reactant
		auto productReactant = network->getCompound(typeName,
				productCompositionVector);
		// If the product exists, mark the proper reaction arrays and add it to the list
		if (productReactant) {
			// Setup the connectivity array for the second reactant
			secondId = secondReactant->getId();
			setReactionConnectivity(secondId);
			// Creates the combining cluster
			// The reaction constant will be computed later and is set to 0.0 for now
			CombiningCluster combCluster(secondReactant, 0.0);
			// Push the product onto the list of clusters that combine with this one
			combiningReactants.push_back(combCluster);
		}
	}

	return;
}

void HeInterstitialCluster::createReactionConnectivity() {
	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);

	// This cluster is always (He_a)(I_b)

	// Helium absorption by HeI clusters producing this cluster
	// He_(a-i) + (He_i)(I_b) --> (He_a)(I_b)
	// Get all the He clusters from the network
	auto reactants = network->getAll(heType);
	auto reactantsSize = reactants.size();
	for (int i = 0; i < reactantsSize; i++) {
		auto heliumReactant = (PSICluster *) reactants[i];
		auto heliumReactantSize = heliumReactant->getSize();
		// Get the second reactant, i.e. HeI cluster with He number smaller
		// by the size of the helium reactant
		auto comp = getComposition();
		std::vector<int> compositionVec = { comp[heType] - heliumReactantSize,
				comp[vType], comp[iType] };
		auto secondReactant = (PSICluster *) network->getCompound(typeName, compositionVec);
		// Create a ReactingPair with the two reactants if they both exist
		if (secondReactant) {
			// The reaction constant will be computed later, it is set to 0.0 for now
			ClusterPair pair(heliumReactant, secondReactant, 0.0);
			// Add the pair to the list
			reactingPairs.push_back(pair);
			// Setup the connectivity array
			int Id = heliumReactant->getId();
			setReactionConnectivity(Id);
			Id = secondReactant->getId();
			setReactionConnectivity(Id);
		}
	}

	// Single Interstitial absorption by HeI clusters producing this cluster
	// (He_a)[I_(b-1)] + I --> (He_a)(I_b)
	// Get the single Interstitial cluster
	auto singleIReactant = (PSICluster *) network->get(iType, 1);
	// Get the second reactant, i.e. HeI cluster with one less I
	auto comp = getComposition();
	std::vector<int> compositionVec = { comp[heType], comp[vType],
			comp[iType] - 1 };
	auto secondReactant = (PSICluster *) network->getCompound(typeName, compositionVec);
	// Create a ReactingPair with the two reactants if they both exist
	if (singleIReactant && secondReactant) {
		// The reaction constant will be computed later, it is set to 0.0 for now
		ClusterPair pair(singleIReactant, secondReactant, 0.0);
		// Add the pair to the list
		reactingPairs.push_back(pair);
		// Setup the connectivity array
		int Id = singleIReactant->getId();
		setReactionConnectivity(Id);
		Id = secondReactant->getId();
		setReactionConnectivity(Id);
	}

	// Helium-Interstitial clustering producing this cluster
	// He_a + I_b --> (He_a)(I_b)
	// Get the He reactant
	auto heliumReactant = (PSICluster *) network->get(heType, numHe);
	// Get the I reactant
	auto interstitialReactant = (PSICluster *) network->get(iType, numI);
	// Create a ReactingPair with the two reactants if they both exist
	if (heliumReactant && interstitialReactant) {
		// The reaction constant will be computed later, it is set to 0.0 for now
		ClusterPair pair(heliumReactant, interstitialReactant, 0.0);
		// Add the pair to the list
		reactingPairs.push_back(pair);
		// Setup the connectivity array
		int Id = heliumReactant->getId();
		setReactionConnectivity(Id);
		Id = interstitialReactant->getId();
		setReactionConnectivity(Id);
	}

	// Interstitial reduction by Vacancy absorption in HeI clusters producing this cluster
	// (He_a)[I_(b+i)] + (V_i) --> (He_a)(I_b)
	// Get all the V clusters from the network
	reactants = network->getAll(vType);
	reactantsSize = reactants.size();
	for (int i = 0; i < reactantsSize; i++) {
		auto vacancyReactant = (PSICluster *) reactants[i];
		auto vacancyReactantSize = vacancyReactant->getSize();
		// Get the second reactant, i.e. HeI cluster with I number bigger
		// by the size of the vacancy reactant
		auto comp = getComposition();
		std::vector<int> compositionVec = { comp[heType],
				comp[vType], comp[iType] + vacancyReactantSize};
		auto secondReactant = (PSICluster *) network->getCompound(typeName, compositionVec);
		// Create a ReactingPair with the two reactants if they both exist
		if (secondReactant) {
			// The reaction constant will be computed later, it is set to 0.0 for now
			ClusterPair pair(vacancyReactant, secondReactant, 0.0);
			// Add the pair to the list
			reactingPairs.push_back(pair);
			// Setup the connectivity array
			int Id = vacancyReactant->getId();
			setReactionConnectivity(Id);
			Id = secondReactant->getId();
			setReactionConnectivity(Id);
		}
	}

	// Interstitial reduction by Vacancy absorption in HeI clusters
	// (He_a)(I_b) + (V_c) --> (He_a)[I_(b-c)]
	// Get all the V clusters from the network
	reactants = network->getAll(vType);
	// replaceInCompound handles this reaction, it is overridden in this class
	replaceInCompound(reactants, iType, vType);

	// Helium absorption by HeI clusters
	// He_c + (He_a)(I_b) --> [He_(a+c)](I_b)
	// Get all the He clusters from the network
	reactants = network->getAll(heType);
	// combineClusters handles He combining with HeI to form HeI
	combineClusters(reactants, typeName);

	// Single Interstitial absorption by HeI clusters
	// (He_a)(I_b) + I --> (He_a)[I_(b+1)]
	// The single Interstitial cluster is already set
	if (singleIReactant) {
		// Create a container for it
		std::vector<Reactant *> singleIInVector;
		singleIInVector.push_back(singleIReactant);
		// Call the combination function even though there is only one cluster
		// because it handles all of the work to properly connect the three
		// clusters in the reaction.
		combineClusters(singleIInVector, typeName);
	}

	return;
}

void HeInterstitialCluster::createDissociationConnectivity() {

	// This cluster is always (He_a)(I_b)

	// He Dissociation
	// (He_a)(I_b) --> [He_(a-1)](I_b) + He
	// Get the cluster with one less helium
	std::vector<int> compositionVec = { numHe - 1, 0, numI };
	auto heIClusterLessHe = (PSICluster *) network->getCompound(typeName, compositionVec);
	// Get the single helium cluster
	auto singleCluster = (PSICluster *) network->get(heType, 1);
	emitClusters(singleCluster, heIClusterLessHe);
	// [He_(a+1)](V_b) --> (He_a)(V_b) + He
	// Get the cluster with one more helium
	compositionVec = { numHe + 1, 0, numI };
	auto heIClusterMoreHe = (PSICluster *) network->getCompound(typeName, compositionVec);
	// Here it is important that heVClusterMoreHe is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heIClusterMoreHe, singleCluster);

	// Interstitial Dissociation
	// (He_a)(I_b) --> He_(a)[I_(b-1)] + I
	// Get the cluster with one less interstitial
	compositionVec = {numHe, 0, numI - 1};
	auto heIClusterLessI = (PSICluster *) network->getCompound(typeName, compositionVec);
	//Get the single interstitial cluster
	singleCluster = (PSICluster *) network->get(iType, 1);
	emitClusters(singleCluster, heIClusterLessI);
	// He_(a)[I_(b+1)] --> (He_a)(I_b) + I
	// Get the cluster with one more interstitial
	compositionVec = {numHe, 0, numI + 1};
	auto heIClusterMoreI = (PSICluster *) network->getCompound(typeName, compositionVec);
	// Here it is important that heIClusterMoreI is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heIClusterMoreI, singleCluster);

	return;
}
