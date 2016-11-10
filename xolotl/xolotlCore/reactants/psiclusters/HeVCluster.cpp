// Includes
#include "HeVCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>

using namespace xolotlCore;

HeVCluster::HeVCluster(int numHe, int numV,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(registry), numHe(numHe), numV(numV) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = numHe + numV;

	// Update the composition map
	compositionMap[heType] = numHe;
	compositionMap[vType] = numV;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << numHe << "V_" << numV;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = heVType;

	// Compute the reaction radius
	reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::tungstenLatticeConstant
			+ pow(
					(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0) * numV)
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
			- pow(
					(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0))
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0));

	return;
}

HeVCluster::HeVCluster(HeVCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numV = other.numV;

	return;
}

void HeVCluster::replaceInCompound(std::vector<IReactant *> & reactants,
		const std::string& oldComponentName) {
	// Local Declarations
	std::map<std::string, int> myComp = getComposition(), productReactantComp;
	int myComponentNumber = myComp[oldComponentName];
	int secondId = 0;

	// Loop over all of the extra reactants in this reaction and handle the replacement
	for (unsigned int i = 0; i < reactants.size(); i++) {
		// Get the second reactant and its size
		auto secondReactant = (PSICluster *) reactants[i];
		int secondReactantSize = secondReactant->getSize();
		// Create the composition vector
		productReactantComp = myComp;
		// Updated the modified components
		productReactantComp[oldComponentName] = myComponentNumber
				- secondReactantSize;
		// Create the composition vector -- FIXME! This should be general!
		std::vector<int> productCompositionVector = {
				productReactantComp[heType], productReactantComp[vType], 0 };
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

void HeVCluster::combineClusters(std::vector<IReactant *> & clusters,
		const std::string& productName) {
	// Initial declarations
	std::map<std::string, int> myComposition = getComposition(),
			secondComposition;

	// Loop on the potential combining reactants
	for (unsigned int i = 0; i < clusters.size(); i++) {
		// Get the second reactant, its composition and its index
		auto secondCluster = (PSICluster *) clusters[i];
		secondComposition = secondCluster->getComposition();
		// Check that the simple product [He_(a+c)](V_b) doesn't exist
		std::vector<int> comp = { myComposition[heType]
				+ secondComposition[heType], myComposition[vType]
				+ secondComposition[vType], myComposition[iType]
				+ secondComposition[iType] };
		auto simpleProduct = network->getCompound(productName, comp);
		if (simpleProduct)
			continue;
		// The simple product doesn't exist so it will go though trap-mutation
		// The reaction is
		// (He_a)(V_b) + He_c --> [He_(a+c)][V_(b+1)] + I
		comp = {myComposition[heType] + secondComposition[heType],
			myComposition[vType] + secondComposition[vType] + 1,
			myComposition[iType] + secondComposition[iType]};
		auto firstProduct = network->getCompound(productName, comp);
		auto secondProduct = network->get(iType, 1);
		// If both products exist
		if (firstProduct && secondProduct) {
			// This cluster combines with the second reactant
			setReactionConnectivity(secondCluster->getId());
			// Creates the combining cluster
			// The reaction constant will be computed later and is set to 0.0 for now
			CombiningCluster combCluster(secondCluster, 0.0);
			// Push the product onto the list of clusters that combine with this one
			combiningReactants.push_back(combCluster);
		}

		// Case with I_2
		// (He_a)(V_b) + He_c --> [He_(a+c)][V_(b+2)] + I_2
		// If [He_(a+c)][V_(b+1)] does not exist
		if (firstProduct)
			continue;

		// Get the new products [He_(a+c)][V_(b+2)] and I_2
		comp = {myComposition[heType] + secondComposition[heType],
			myComposition[vType] + secondComposition[vType] + 2,
			myComposition[iType] + secondComposition[iType]};
		firstProduct = network->getCompound(productName, comp);
		secondProduct = network->get(iType, 2);
		// If both products exist
		if (firstProduct && secondProduct) {
			// This cluster combines with the second reactant
			setReactionConnectivity(secondCluster->getId());
			// Creates the combining cluster
			// The reaction constant will be computed later and is set to 0.0 for now
			CombiningCluster combCluster(secondCluster, 0.0);
			// Push the product onto the list of clusters that combine with this one
			combiningReactants.push_back(combCluster);
		}
	}

	return;
}

void HeVCluster::createReactionConnectivity() {
	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);

	// This cluster is always (He_a)(V_b)

	// Helium absorption by HeV clusters producing this cluster
	// He_(a-i) + (He_i)(V_b) --> (He_a)(V_b)
	// Get all the He clusters from the network
	auto reactants = network->getAll(heType);
	for (unsigned int i = 0; i < reactants.size(); i++) {
		auto heliumReactant = (PSICluster *) reactants[i];
		int heliumReactantSize = heliumReactant->getSize();
		// Get the second reactant, i.e. HeV cluster with He number smaller
		// by the size of the helium reactant
		auto comp = getComposition();
		std::vector<int> compositionVec = { comp[heType] - heliumReactantSize,
				comp[vType], 0 };
		auto secondReactant = (PSICluster *) network->getCompound(typeName,
				compositionVec);
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

	// Single Vacancy absorption by HeV clusters producing this cluster
	// (He_a)[V_(b-1)] + V --> (He_a)(V_b)
	// Get the single Vacancy cluster
	auto singleVReactant = (PSICluster *) network->get(vType, 1);
	// Get the second reactant, i.e. HeV cluster with one less V
	auto comp = getComposition();
	std::vector<int> compositionVec = { comp[heType], comp[vType] - 1, 0 };
	auto secondReactant = (PSICluster *) network->getCompound(typeName,
			compositionVec);
	// Create a ReactingPair with the two reactants if they both exist
	if (singleVReactant && secondReactant) {
		// The reaction constant will be computed later, it is set to 0.0 for now
		ClusterPair pair(singleVReactant, secondReactant, 0.0);
		// Add the pair to the list
		reactingPairs.push_back(pair);
		// Setup the connectivity array
		int Id = singleVReactant->getId();
		setReactionConnectivity(Id);
		Id = secondReactant->getId();
		setReactionConnectivity(Id);
	}

	// Helium-Vacancy clustering producing this cluster
	// He_a + V_b --> (He_a)(V_b)
	// Get the He reactant
	auto heliumReactant = (PSICluster *) network->get(heType, numHe);
	// Get the V reactant
	auto vacancyReactant = (PSICluster *) network->get(vType, numV);
	// Create a ReactingPair with the two reactants if they both exist
	if (heliumReactant && vacancyReactant) {
		// The reaction constant will be computed later, it is set to 0.0 for now
		ClusterPair pair(heliumReactant, vacancyReactant, 0.0);
		// Add the pair to the list
		reactingPairs.push_back(pair);
		// Setup the connectivity array
		int Id = heliumReactant->getId();
		setReactionConnectivity(Id);
		Id = vacancyReactant->getId();
		setReactionConnectivity(Id);
	}

	// Vacancy reduction by Interstitial absorption in HeV producing this cluster
	// (He_a)[V_(b+c)] + (I_c) --> (He_a)(V_b)
	// Get all the I clusters from the network
	reactants = network->getAll(iType);
	// Get the composition of this cluster
	comp = getComposition();
	for (unsigned int i = 0; i < reactants.size(); i++) {
		auto interstitialReactant = (PSICluster *) reactants[i];
		int interstitialReactantSize = interstitialReactant->getSize();
		// Get the second reactant, i.e. HeV cluster with V number bigger
		// by the size of the interstitial reactant
		std::vector<int> compositionVec = { comp[heType], comp[vType]
				+ interstitialReactantSize, 0 };
		auto secondReactant = (PSICluster *) network->getCompound(typeName,
				compositionVec);
		// Create a ReactingPair with the two reactants if they both exist
		if (secondReactant) {
			// The reaction constant will be computed later, it is set to 0.0 for now
			ClusterPair pair(interstitialReactant, secondReactant, 0.0);
			// Add the pair to the list
			reactingPairs.push_back(pair);
			// Setup the connectivity array
			int Id = interstitialReactant->getId();
			setReactionConnectivity(Id);
			Id = secondReactant->getId();
			setReactionConnectivity(Id);
		}
	}

	// Vacancy reduction by Interstitial absorption in HeV
	// (He_a)(V_b) + (I_c) --> (He_a)[V_(b-c)]
	// Get all the I clusters from the network
	reactants = network->getAll(iType);
	// replaceInCompound handles this reaction, it is overridden in this class
	replaceInCompound(reactants, vType);

	// Helium absorption by HeV clusters
	// He_c + (He_a)(V_b) --> [He_(a+c)](V_b)
	// Get all the He clusters from the network
	reactants = network->getAll(heType);
	// combineClusters handles He combining with HeV to form HeV
	PSICluster::combineClusters(reactants, typeName);

	// Single Vacancy absorption by HeV clusters
	// (He_a)(V_b) + V --> (He_a)[V_(b+1)]
	// The single Vacancy cluster is already set
	if (singleVReactant) {
		// Create a container for it
		std::vector<IReactant *> singleVInVector;
		singleVInVector.push_back(singleVReactant);
		// Call the combination function even though there is only one cluster
		// because it handles all of the work to properly connect the three
		// clusters in the reaction.
		PSICluster::combineClusters(singleVInVector, typeName);
	}

	// Helium absorption leading to trap mutation
	// (He_a)(V_b) + He_c --> [He_(a+c)][V_(b+1)] + I
	// or
	// (He_a)(V_b) + He_c --> [He_(a+c)][V_(b+2)] + I_2
	// HeVCluster::combineClusters handles He combining with HeV to go through trap-mutation
	combineClusters(reactants, typeName);

	// Helium absorption by HeV cluster leading to trap mutation and the production of this cluster
	// [He_(a-c)][V_(b-1)] + He_c --> (He_a)(V_b) + I
	// or
	// [He_(a-c)][V_(b-2)] + He_c --> (He_a)(V_b) + I_2
	// Happens only if (He_a)[V_(b-1)] is not present in the network
	PSICluster * smallerCluster;
	// (b-1) can be 0 so (He_a)[V_(b-1)] can be a helium cluster
	if (comp[vType] == 1) {
		smallerCluster = (PSICluster *) network->get(heType, comp[heType]);
	} else {
		std::vector<int> compositionVec = { comp[heType], comp[vType] - 1, 0 };
		smallerCluster = (PSICluster *) network->getCompound(typeName,
				compositionVec);
	}
	if (!smallerCluster) {
		// Loop on the possible He reactants He_c
		for (unsigned int i = 0; i < reactants.size(); i++) {
			auto heReactant = (PSICluster *) reactants[i];

			// Case with I
			// Get the other reactant [He_(a-c)][V_(b-1)] that can be He or HeV
			PSICluster * otherReactant;
			if (comp[vType] == 1) {
				// We want (a-c) to be smaller or equal to c in order to avoid double counting
				if (comp[heType] > 2 * heReactant->getSize())
					continue;

				otherReactant = (PSICluster *) network->get(heType,
						comp[heType] - heReactant->getSize());
			} else {
				std::vector<int> compositionVec = { comp[heType]
						- heReactant->getSize(), comp[vType] - 1, 0 };
				otherReactant = (PSICluster *) network->getCompound(typeName,
						compositionVec);
			}
			// Get I the other product
			auto iCluster = (PSICluster *) network->get(iType, 1);
			// If the other reactant and product exist
			if (otherReactant && iCluster) {
				// The reaction is really allowed
				// Create the pair
				// The reaction constant will be computed later, it is set to 0.0 for now
				ClusterPair pair(heReactant, otherReactant, 0.0);
				// Add the pair to the list
				reactingPairs.push_back(pair);
				// Setup the connectivity array
				int Id = heReactant->getId();
				setReactionConnectivity(Id);
				Id = otherReactant->getId();
				setReactionConnectivity(Id);
			}

			// Case with I_2
			// Get the other reactant [He_(a-c)][V_(b-2)] that can be He or HeV
			if (comp[vType] == 2) {
				// We want (a-c) to be smaller or equal to c in order to avoid double counting
				if (comp[heType] > 2 * heReactant->getSize())
					continue;

				otherReactant = (PSICluster *) network->get(heType,
						comp[heType] - heReactant->getSize());
			} else {
				std::vector<int> compositionVec = { comp[heType]
						- heReactant->getSize(), comp[vType] - 2, 0 };
				otherReactant = (PSICluster *) network->getCompound(typeName,
						compositionVec);
			}
			// Get I_2 the other product
			iCluster = (PSICluster *) network->get(iType, 2);
			// If the other reactant and product exist
			if (otherReactant && iCluster) {
				// The reaction is really allowed
				// Create the pair
				// The reaction constant will be computed later, it is set to 0.0 for now
				ClusterPair pair(heReactant, otherReactant, 0.0);
				// Add the pair to the list
				reactingPairs.push_back(pair);
				// Setup the connectivity array
				int Id = heReactant->getId();
				setReactionConnectivity(Id);
				Id = otherReactant->getId();
				setReactionConnectivity(Id);
			}
		}
	}

	return;
}

void HeVCluster::createDissociationConnectivity() {
	// This cluster is always (He_a)(V_b)

	// He Dissociation
	// (He_a)(V_b) --> [He_(a-1)](V_b) + He
	// Get the cluster with one less helium
	std::vector<int> compositionVec = { numHe - 1, numV, 0 };
	auto heVClusterLessHe = (PSICluster *) network->getCompound(typeName,
			compositionVec);
	// Special case for numHe = 1
	if (numHe == 1) {
		heVClusterLessHe = (PSICluster *) network->get(vType, numV);
	}
	// Get the single helium cluster
	auto singleCluster = (PSICluster *) network->get(heType, 1);
	emitClusters(singleCluster, heVClusterLessHe);
	// [He_(a+1)](V_b) --> (He_a)(V_b) + He
	// Get the cluster with one more helium
	compositionVec = {numHe + 1, numV, 0};
	auto heVClusterMoreHe = (PSICluster *) network->getCompound(typeName,
			compositionVec);
	// Here it is important that heVClusterMoreHe is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heVClusterMoreHe, singleCluster);

	// Vacancy Dissociation
	// (He_a)(V_b) --> He_(a)[V_(b-1)] + V
	// Get the cluster with one less vacancy
	compositionVec = {numHe, numV - 1, 0};
	auto heVClusterLessV = (PSICluster *) network->getCompound(typeName,
			compositionVec);
	// Special case for numV = 1
	if (numV == 1) {
		heVClusterLessV = (PSICluster *) network->get(heType, numHe);
	}
	// Skip He_1V_1 because it was counted in the He dissociation
	if (numHe == 1 && numV == 1)
		heVClusterLessV = nullptr;

	// Get the single vacancy cluster
	singleCluster = (PSICluster *) network->get(vType, 1);
	emitClusters(singleCluster, heVClusterLessV);
	// He_(a)[V_(b+1)] --> (He_a)(V_b) + V
	// Get the cluster with one more vacancy
	compositionVec = {numHe, numV + 1, 0};
	auto heVClusterMoreV = (PSICluster *) network->getCompound(typeName,
			compositionVec);
	// Here it is important that heVClusterMoreV is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heVClusterMoreV, singleCluster);

	return;
}
