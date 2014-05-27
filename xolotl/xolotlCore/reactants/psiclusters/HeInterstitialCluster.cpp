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

}

HeInterstitialCluster::HeInterstitialCluster(const HeInterstitialCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numI = other.numI;
}

HeInterstitialCluster::~HeInterstitialCluster() {
}

std::shared_ptr<Reactant> HeInterstitialCluster::clone() {
	std::shared_ptr<Reactant> reactant = std::make_shared<HeInterstitialCluster>(*this);
	return reactant;
}

double HeInterstitialCluster::getGenByEm() {
	return 0;
}

double HeInterstitialCluster::getAnnByEm() {
	return 0;
}

void HeInterstitialCluster::createReactionConnectivity() {

	// Local Declarations
	auto psiNetwork = std::dynamic_pointer_cast<PSIClusterReactionNetwork>(
			network);
	auto props = psiNetwork->getProperties();
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	int maxHeIClusterSize = std::stoi(props["maxHeIClusterSize"]);
	std::vector<int> firstComposition, secondComposition;

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(getId());

	/* ----- (He_a)(I_b) + (V_c) --> (He_a)[I_(b-c)] -----
	 * This section adds the clusters that produce this cluster to the array
	 * for vacancy absorption by HeI.
	 */
	for (int c = 1; c <= maxVClusterSize; c++) { // hehe... c++!
		// Get the first reactant's composition and then retrieve it
		firstComposition = psiNetwork->getCompositionVector(numHe, 0, numI + c);
		auto firstReactant = psiNetwork->getCompound("HeI", firstComposition);
		// Set the second reactant
		auto secondReactant = psiNetwork->get(vType, c);
		// Create a ReactingPair with the two reactants
		if (firstReactant && secondReactant) {
			ReactingPair pair((PSICluster *) firstReactant, (PSICluster *) secondReactant);
			// Add the pair to the list
			reactingPairs.push_back(pair);
		}
	}

	/* ----- (He_a)(I_b) + I --> (He_a)[I_(b+1)] -----
	 * This section adds the clusters that produce this cluster to the array
	 * for single species interstitial absorption by HeI.
	 *
	 * This section also handles the case of this cluster combining with a
	 * single interstitial to produce an HeI cluster with one addition
	 * interstitial.
	 */
	firstComposition = psiNetwork->getCompositionVector(numHe, 0, numI - 1);
	// Get those Reactants from the network
	auto firstReactant = psiNetwork->getCompound("HeI", firstComposition);
	auto secondReactant = psiNetwork->get(iType, 1);
	if (firstReactant && secondReactant) {
		// Create the Reacting Pair
		ReactingPair pair((PSICluster *) firstReactant, (PSICluster *) secondReactant);
		// Add the pair to the list
		reactingPairs.push_back(pair);
		// Add single I to the list of clusters this one interacts with if it
		// doesn't violate the maximum size limit.
		setReactionConnectivity(secondReactant->getId());
		combiningReactants.push_back(secondReactant);
	}

	/* ----- (He_a)(I_b) + He_c --> [He_(a+c)](I_b)
	 * HeI clusters can absorb helium clusters of any size so long as the
	 * maximum size limit is not violated.
	 */
	auto reactants = psiNetwork->getAll(heType);
	combineClusters(reactants, maxHeIClusterSize, "HeI");

	/* ----- (A*He)(B*I) + (C*V) --> (A*He)[(B-C)*I] -----
	 * This section adds the clusters that are produced by this cluster to the
	 * array for vacancy absorption by HeI.
	 */
	reactants = psiNetwork->getAll(vType);
	replaceInCompound(reactants, vType, iType);

	// Set the references to the size one clusters
	heCluster = (PSICluster *) network->get(heType, 1);
	vCluster = (PSICluster *) network->get(vType, 1);
	iCluster = (PSICluster *) network->get(iType, 1);

	// Set the references for clusters one size smaller than this cluster,
	// starting with one less He.
	std::vector<int> compositionVec = { numHe - 1, 0, numI };
	heIClusterLessHe = (PSICluster *) network->getCompound(heIType, compositionVec);
	// Store the cluster with one less vacancy
	compositionVec = {numHe, 0, numI - 1};
	heIClusterLessI = (PSICluster *) network->getCompound(heIType, compositionVec);

	return;
}

void HeInterstitialCluster::createDissociationConnectivity() {

	// Local Declarations
	auto psiNetwork = std::dynamic_pointer_cast<PSIClusterReactionNetwork>(
			network);
	auto props = psiNetwork->getProperties();
	std::vector<int> composition;

	// Get the required dissociating clusters. These are stored for the flux
	// computation later.
	auto heCluster = network->get(heType, 1);
	auto vCluster = network->get(vType, 1);
	auto iCluster = network->get(iType, 1);

	// Store the cluster with one less helium
	std::vector<int> compositionVec = { numHe - 1, 0, numI };
	auto heIClusterLessHe = network->getCompound("HeI", compositionVec);
	// Store the cluster with one less vacancy
	compositionVec = {numHe, 0, numI - 1};
	auto heIClusterLessI = network->getCompound("HeI", compositionVec);

	// He Dissociation, get the [(numHe-1)*He]I and He
	composition = psiNetwork->getCompositionVector(numHe - 1, 0, numI);
	auto otherMixedCluster = psiNetwork->getCompound("HeI", composition);
	auto singleCluster = psiNetwork->get(heType, 1);
	dissociateClusters(singleCluster, otherMixedCluster);

	// Interstitial Dissociation, get He[(numV-1)*I] and I
	composition = psiNetwork->getCompositionVector(numHe, 0, numI - 1);
	otherMixedCluster = psiNetwork->getCompound("HeI", composition);
	singleCluster = psiNetwork->get(iType, 1);
	dissociateClusters(singleCluster, otherMixedCluster);

	return;
}

void HeInterstitialCluster::setTemperature(double temp) {

	// Call the base class version to set all of the basic quantities.
	PSICluster::setTemperature(temp);

	// Calculate the much easier f4 term... first
	f4 = calculateDissociationConstant(*this, *heCluster, temperature)
			+ calculateDissociationConstant(*this, *vCluster, temperature)
			+ calculateDissociationConstant(*this, *iCluster, temperature);

	return;
}

double HeInterstitialCluster::getDissociationFlux(double temperature) const {

	// Local Declarations
	double f3 = 0.0;

	// Only dissociate if possible
	if (heCluster && vCluster && iCluster) {
		// FIXME! Make sure that this works as expected! Make sure that it
		// correctly picks out right component in
		// calculateDissociationConstant!

		// Loop over all the elements of the dissociation
		// connectivity to find where this mixed species dissociates and
		// calculate the f3 term.
		//
		// TODO - What's the performance difference between getting all of the
		// reactants and pulling each reactant separately?
		auto reactants = network->getAll();
		auto dissociatingSet = getDissociationConnectivitySet();
		for (auto it = dissociatingSet.begin(); it != dissociatingSet.end();
				it++) {
			// Set the current reactant
			auto currentCluster = (PSICluster *) reactants->at(*it-1);
			// Get the cluster map of this connection
			auto dissClusterComposition = currentCluster->getComposition();
			// We need to find if this is a Helium dissociation
			PSICluster * secondCluster;
			if (numHe - dissClusterComposition[heType] == 1 && numI == dissClusterComposition[iType]
					&& dissClusterComposition[vType] == 0) {
				secondCluster = heCluster;
			} else if (numHe == dissClusterComposition[heType]
					&& numI - dissClusterComposition[iType] == 1 && dissClusterComposition[vType] == 0) {
				// trap mutation
				secondCluster = vCluster;
			} else if (numHe == dissClusterComposition[heType]
					&& dissClusterComposition[vType] - numI == 1 && dissClusterComposition[iType] == 0) {
				// or interstitial dissociation
				secondCluster = iCluster;
			}
			// Update the flux calculation
			if (secondCluster) {
				// Get the single species cluster that comes out with it
				f3 += calculateDissociationConstant(*currentCluster,
						*secondCluster, temperature)
						* currentCluster->getConcentration();
			}
		}
	}

	return f3 - (f4 * getConcentration());
}

/**
 * This operation computes the partial derivatives due to dissociation
 * reactions. The partial derivatives due to dissociation for compound
 * clusters are significantly different than those single-species clusters.
 *
 * @param partials The vector into which the partial derivatives should be
 * inserted. This vector should have a length equal to the size of the
 * network.
 * @param temperature The temperature at which the reactions are occurring.
 */
void HeInterstitialCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials, double temperature) const {

	// Local Declarations
	int index = 0;

	// Partial derivative with respect to changes in this cluster
	double partialDeriv = calculateDissociationConstant(*this, *heCluster,
			temperature)
			+ calculateDissociationConstant(*this, *vCluster, temperature)
			+ calculateDissociationConstant(*this, *iCluster, temperature);
	// Add it to the list of partials
	partials[index] += partialDeriv;

	// Compute the partial derivative if the smaller clusters exists
	if (heIClusterLessHe) {
		partialDeriv = calculateDissociationConstant(*heIClusterLessHe,
				*(PSICluster *)heCluster, temperature);
		index = heIClusterLessHe->getId() - 1;
		partials[index] += partialDeriv;
	}

	// Compute the partial derivative if the smaller clusters exists
	if (heIClusterLessI) {
		partialDeriv = calculateDissociationConstant(*heIClusterLessI,
				*(PSICluster *)vCluster, temperature);
		index = heIClusterLessI->getId() - 1;
		partials[index] += partialDeriv;
	}

	// This cluster cannot dissociate into a smaller HeI cluster and a vacancy,
	// so there is no partial derivative term for that case.

	return;
}
