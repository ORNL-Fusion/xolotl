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

void HeVCluster::createReactionConnectivity() {

	// Local Declarations
	auto psiNetwork = std::dynamic_pointer_cast<PSIClusterReactionNetwork>(
			network);
	auto props = psiNetwork->getProperties();
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxHeVClusterSize = std::stoi(props["maxHeVClusterSize"]);
	std::vector<int> firstComposition, secondComposition;

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(getId());

	/* ----- (He_a)(V_b) + (He_c) --> [He_(a+c)]*(V_b) -----
	 * Helium absorption by HeV clusters that results
	 * in the production of this cluster.
	 */
	for (int z = 1; z <= maxHeClusterSize; z++) {
		// Get the first reactant
		firstComposition = psiNetwork->getCompositionVector(numHe - z, numV, 0);
		auto firstReactant = psiNetwork->getCompound("HeV", firstComposition);
		// Get the second reactant
		auto secondReactant = psiNetwork->get(heType, z);
		// Create a ReactingPair with the two reactants
		if (firstReactant && secondReactant) {
			ReactingPair pair((PSICluster *) firstReactant, (PSICluster *) secondReactant);
			// Add the pair to the list
			reactingPairs.push_back(pair);
		}
	}

	/* ----- (He_a)(V_b) + V --> (He_a)[V_(b+1)] -----
	 * HeV clusters are also produced by single-vacancy absorption by another
	 * HeV cluster. In this case, (A*He)[(B-1)*V] produces the current cluster.
	 */
	firstComposition = psiNetwork->getCompositionVector(numHe, numV - 1, 0);
	auto firstReactant = psiNetwork->getCompound("HeV", firstComposition);
	auto secondReactant = psiNetwork->get(vType, 1);
	// Create a ReactingPair with the two reactants
	if (firstReactant && secondReactant) {
		ReactingPair pair((PSICluster *) firstReactant, (PSICluster *) secondReactant);
		// Add the pair to the list
		reactingPairs.push_back(pair);
	}

	/* ----- (He_a) + (V_b) --> (He_a)(V_b) -----
	 * Helium-vacancy clustering that results
	 * in the production of this cluster.
	 */
	// Get the first reactant
	firstReactant = psiNetwork->get(heType, numHe);
	// Get the second reactant
	secondReactant = psiNetwork->get(vType, numV);
	// Create a ReactingPair with the two reactants
	if (firstReactant && secondReactant) {
		ReactingPair pair((PSICluster *) firstReactant, (PSICluster *) secondReactant);;
		// Add the pair to the list
		reactingPairs.push_back(pair);
	}

	/* ----- (He_a)*(V_b) + I_c  --> (He_a)*[V_(b-c)] -----
	 * Helium-vacancy clusters lose vacancies when they interact with
	 * interstitial clusters.
	 *
	 * We assume that the HeV and interstitial cluster can only
	 * interact if they would produce another HeV cluster, not single He.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the (outgoing) flux due to combination
	 * reactions.
	 */
	auto reactants = psiNetwork->getAll(iType);
	replaceInCompound(reactants, iType, vType);

	/* ---- (He_a)*(V_b) + He_c --> [He_(a+c)]*(V_b) ----
	 * HeV clusters can absorb helium clusters so long as they do not cross
	 * the max size limit.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	reactants = psiNetwork->getAll(heType);
	combineClusters(reactants, maxHeVClusterSize, "HeV");

	/* ----- (He_a)*(V_b) + V --> (He_a)*[V_(b+1)] -----
	 * HeV clusters can absorb single vacancies.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	secondReactant = psiNetwork->get(vType, 1);
	if (secondReactant) {
		// Create a container for it
		std::vector<Reactant *> singleVInVector;
		singleVInVector.push_back(secondReactant);
		// Call the combination function even though there is only one cluster
		// because it handles all of the work to properly connect the three
		// clusters in the reaction.
		combineClusters(singleVInVector,maxHeVClusterSize,"HeV");
	}

	// Set the references to the size one clusters
	heCluster = (PSICluster *) network->get(heType, 1);
	vCluster = (PSICluster *) network->get(vType, 1);
	iCluster = (PSICluster *) network->get(iType, 1);

	// Set the references for clusters one size smaller than this cluster,
	// starting with one less He.
	std::vector<int> compositionVec = { numHe - 1, numV, 0 };
	heVClusterLessHe = (PSICluster *) network->getCompound(heVType, compositionVec);
	// Store the cluster with one less vacancy
	compositionVec = {numHe, numV - 1, 0};
	heVClusterLessV = (PSICluster *) network->getCompound(heVType, compositionVec);

	return;
}

void HeVCluster::createDissociationConnectivity() {

	// Local Declarations
	auto psiNetwork = std::dynamic_pointer_cast<PSIClusterReactionNetwork>(
			network);
	auto props = psiNetwork->getProperties();
	std::vector<int> composition;

	// Store the cluster with one less helium
	std::vector<int> compositionVec = { numHe - 1, numV, 0 };
	auto heVClusterLessHe = network->getCompound("HeV", compositionVec);
	// Store the cluster with one less vacancy
	compositionVec = {numHe, numV - 1, 0};
	auto heVClusterLessV = network->getCompound("HeV", compositionVec);

	// He Dissociation, get the [(numHe-1)*He]V and He
	composition = psiNetwork->getCompositionVector(numHe - 1, numV, 0);
	auto otherMixedCluster = psiNetwork->getCompound("HeV", composition);
	auto singleCluster = psiNetwork->get(heType, 1);
	dissociateClusters(singleCluster, otherMixedCluster);

	// Vacancy Dissociation, get He[(numV-1)*V] and V
	composition = psiNetwork->getCompositionVector(numHe, numV - 1, 0);
	otherMixedCluster = psiNetwork->getCompound("HeV", composition);
	singleCluster = psiNetwork->get(vType, 1);
	dissociateClusters(singleCluster, otherMixedCluster);

	// Trap mutation, get He[(numV+1)*V] and I
	composition = psiNetwork->getCompositionVector(numHe, numV + 1, 0);
	otherMixedCluster = psiNetwork->getCompound("HeV", composition);
	singleCluster = psiNetwork->get(iType, 1);
	dissociateClusters(singleCluster, otherMixedCluster);

	return;
}

void HeVCluster::setTemperature(double temp) {

	// Call the base class version to set all of the basic quantities.
	PSICluster::setTemperature(temp);

	// Calculate the much easier f4 term... first
	f4 = calculateDissociationConstant(*this, *heCluster, temperature)
			+ calculateDissociationConstant(*this, *vCluster, temperature)
			+ calculateDissociationConstant(*this, *iCluster, temperature);

	return;
}

double HeVCluster::getDissociationFlux(double temperature) const {

	// Local Declarations
	PSICluster * currentCluster;
	double f3 = 0.0;

	// Get the required dissociating clusters. These are stored for the flux
	// computation later.'

	// Only dissociate if possible
	if (heCluster && vCluster && iCluster) {
		// FIXME! Make sure that this works as expected! Make sure that it
		// correctly picks out the right component in
		// calculateDissociationConstant!

		// Loop over all the elements of the dissociation
		// connectivity to find where this mixed species dissociates
		//
		// TODO - What's the performance difference between getting all of the
		// reactants and pulling each reactant separately?
		auto reactants = network->getAll();
		auto dissociatingSet = getDissociationConnectivitySet();
		for (auto it = dissociatingSet.begin(); it != dissociatingSet.end();
				it++) {
			// Set the current reactant
			currentCluster = (PSICluster *) reactants->at(*it - 1);
			// Get the cluster map of this connection
			auto dissClusterComposition = currentCluster->getComposition();
			// We need to find if this is a Helium dissociation
			if (numHe - dissClusterComposition[heType] == 1 && numV == dissClusterComposition[vType]
					&& dissClusterComposition[iType] == 0) {
				f3 += calculateDissociationConstant(*currentCluster,
						*heCluster, temperature)
						* currentCluster->getConcentration();
			} else if (numHe == dissClusterComposition[heType]
					&& numV - dissClusterComposition[vType] == 1 && dissClusterComposition[iType] == 0) {
				// vacancy dissociation
				f3 += calculateDissociationConstant(*currentCluster,
						*vCluster, temperature)
						* currentCluster->getConcentration();
			} else if (numHe == dissClusterComposition[heType]
					&& dissClusterComposition[iType] - numV == 1 && dissClusterComposition[vType] == 0) {
				// or a trap mutation.
				f3 += calculateDissociationConstant(*currentCluster,
						*iCluster, temperature)
						* currentCluster->getConcentration();
			}
		}
	}
	return f3 - f4 * getConcentration();
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
void HeVCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials, double temperature) const {

	// Local Declarations
	int index = 0;

	// Partial derivative with respect to changes in this cluster
	double partialDeriv = calculateDissociationConstant(*this, *heCluster,
			temperature)
			+ calculateDissociationConstant(*this, *vCluster, temperature)
			+ calculateDissociationConstant(*this, *iCluster, temperature);
	// Add it to the list of partials
	partials[thisNetworkIndex] += partialDeriv;

	// Compute the partial derivative if the cluster with one less He exists
	if (heVClusterLessHe) {
		partialDeriv = calculateDissociationConstant(*heVClusterLessHe,
				*heCluster, temperature);
		index = heVClusterLessHe->getId() - 1;
		partials[index] += partialDeriv;
	}

	// Compute the partial derivative if the cluster with one less V exists
	if (heVClusterLessV) {
		partialDeriv = calculateDissociationConstant(*heVClusterLessV,
				*vCluster, temperature);
		index = heVClusterLessV->getId() - 1;
		partials[index] += partialDeriv;
	}

	// This cluster cannot dissociate into a smaller HeV cluster and an
	// interstitial, so there is no partial derivative term for that case.

	return;
}
