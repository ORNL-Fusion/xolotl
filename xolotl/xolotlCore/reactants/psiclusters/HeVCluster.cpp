// Includes
#include "HeVCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>

using namespace xolotlCore;

HeVCluster::HeVCluster(int numHe, int numV,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(1, registry), numHe(numHe), numV(numV) {

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
	typeName = "HeV";

	// Compute the reaction radius
	reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant
			+ pow(
					(3.0 * pow(xolotlCore::latticeConstant, 3.0) * numV)
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
			- pow(
					(3.0 * pow(xolotlCore::latticeConstant, 3.0))
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0));

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

void HeVCluster::replaceInCompound(std::vector<Reactant *> & reactants,
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

void HeVCluster::createReactionConnectivity() {
	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(getId());

	// This cluster is always (He_a)(V_b)

	// Helium absorption by HeV clusters producing this cluster
	// He_(a-i) + (He_i)(V_b) --> (He_a)(V_b)
	// Get all the He clusters from the network
	auto reactants = network->getAll(heType);
	auto reactantsSize = reactants.size();
	for (int i = 0; i < reactantsSize; i++) {
		auto heliumReactant = (PSICluster *) reactants[i];
		auto heliumReactantSize = heliumReactant->getSize();
		// Get the second reactant, i.e. HeV cluster with He number smaller
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

	// Single Vacancy absorption by HeV clusters producing this cluster
	// (He_a)[V_(b-1)] + V --> (He_a)(V_b)
	// Get the single Vacancy cluster
	auto singleVReactant = (PSICluster *) network->get(vType, 1);
	// Get the second reactant, i.e. HeV cluster with one less V
	auto comp = getComposition();
	std::vector<int> compositionVec = { comp[heType], comp[vType] - 1,
			comp[iType] };
	auto secondReactant = (PSICluster *) network->getCompound(typeName, compositionVec);
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
	// (He_a)[V_(b+i)] + (I_i) --> (He_a)(V_b)
	// Get all the I clusters from the network
	reactants = network->getAll(iType);
	reactantsSize = reactants.size();
	for (int i = 0; i < reactantsSize; i++) {
		auto interstitialReactant = (PSICluster *) reactants[i];
		auto interstitialReactantSize = interstitialReactant->getSize();
		// Get the second reactant, i.e. HeV cluster with V number bigger
		// by the size of the interstitial reactant
		auto comp = getComposition();
		std::vector<int> compositionVec = { comp[heType],
				comp[vType] + interstitialReactantSize, comp[iType] };
		auto secondReactant = (PSICluster *) network->getCompound(typeName, compositionVec);
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
	replaceInCompound(reactants, vType, iType);

	// Helium absorption by HeV clusters
	// He_c + (He_a)(V_b) --> [He_(a+c)](V_b)
	// Get all the He clusters from the network
	reactants = network->getAll(heType);
	// combineClusters handles He combining with HeV to form HeV
	combineClusters(reactants, typeName);

	// Single Vacancy absorption by HeV clusters
	// (He_a)(V_b) + V --> (He_a)[V_(b+1)]
	// The single Vacancy cluster is already set
	if (singleVReactant) {
		// Create a container for it
		std::vector<Reactant *> singleVInVector;
		singleVInVector.push_back(singleVReactant);
		// Call the combination function even though there is only one cluster
		// because it handles all of the work to properly connect the three
		// clusters in the reaction.
		combineClusters(singleVInVector, typeName);
	}

	// Set the references to the size one clusters
	heCluster = (PSICluster *) network->get(heType, 1);
	vCluster = (PSICluster *) network->get(vType, 1);
	iCluster = (PSICluster *) network->get(iType, 1);

	return;
}

void HeVCluster::createDissociationConnectivity() {
	// Store the cluster with one less helium
	std::vector<int> compositionVec = { numHe - 1, numV, 0 };
	auto heVClusterLessHe = network->getCompound(typeName, compositionVec);
	// Store the cluster with one more helium
	compositionVec = { numHe + 1, numV, 0 };
	auto heVClusterMoreHe = network->getCompound(typeName, compositionVec);
	// Store the cluster with one less vacancy
	compositionVec = {numHe, numV - 1, 0};
	auto heVClusterLessV = network->getCompound(typeName, compositionVec);
	// Store the cluster with one more vacancy
	compositionVec = {numHe, numV + 1, 0};
	auto heVClusterMoreV = network->getCompound(typeName, compositionVec);

	// This cluster is always (He_a)(V_b)

	// He Dissociation
	// (He_a)(V_b) --> [He_(a-1)](V_b) + He
	auto singleCluster = network->get(heType, 1);
	emitClusters(singleCluster, heVClusterLessHe);
	// [He_(a+1)](V_b) --> (He_a)(V_b) + He
	// Here it is important that heVClusterMoreHe is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heVClusterMoreHe, singleCluster);

	// Vacancy Dissociation
	// (He_a)(V_b) --> He_(a)[V_(b-1)] + V
	singleCluster = network->get(vType, 1);
	emitClusters(singleCluster, heVClusterLessV);
	// He_(a)[V_(b+1)] --> (He_a)(V_b) + V
	// Here it is important that heVClusterMoreV is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heVClusterMoreV, singleCluster);

	// Trap mutation
	// (He_a)(V_b) --> He_(a)[V_(b+1)] + I
	singleCluster = network->get(iType, 1);
	emitClusters(singleCluster, heVClusterMoreV);
	// He_(a)[V_(b-1)] --> (He_a)(V_b) + I
	// Here it is important that heVClusterLessV is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heVClusterLessV, singleCluster);

	// Trap mutation of He
	// He_a --> (He_a)(V_b) + I_b
	// for a = 9 and b = 1
	if (numHe == 9 && numV == 1) {
		// Get He_a
		auto dissociatingReactant = network->get(heType, 9);
		// Get I_b
		auto singleReactant = network->get(iType, 1);
		dissociateCluster(dissociatingReactant, singleReactant);
	}

	return;
}

void HeVCluster::setTemperature(double temp) {

	// Call the base class version to set all of the basic quantities.
	PSICluster::setTemperature(temp);

	return;
}
