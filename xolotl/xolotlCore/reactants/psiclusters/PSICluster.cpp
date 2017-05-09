#include "PSICluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

PSICluster::PSICluster() :
		Reactant() {
	// Set the reactant name appropriately
	name = "PSICluster";
	// Setup the composition map.
	compositionMap[heType] = 0;
	compositionMap[vType] = 0;
	compositionMap[iType] = 0;

	return;
}

PSICluster::PSICluster(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		Reactant(registry) {
	// Set the reactant name appropriately
	name = "PSICluster";
	// Setup the composition map.
	compositionMap[heType] = 0;
	compositionMap[vType] = 0;
	compositionMap[iType] = 0;

	return;
}

// The copy constructor
PSICluster::PSICluster(PSICluster &other) :
		Reactant(other), reactingPairs(other.reactingPairs), combiningReactants(
				other.combiningReactants), dissociatingPairs(
				other.dissociatingPairs), emissionPairs(other.emissionPairs) {
	compositionMap[heType] = other.compositionMap[heType];
	compositionMap[vType] = other.compositionMap[vType];
	compositionMap[iType] = other.compositionMap[iType];

	// Recompute all of the temperature-dependent quantities
	setTemperature(other.getTemperature());

	return;
}

void PSICluster::createProduction(
		std::shared_ptr<ProductionReaction> reaction) {
	// Create a cluster pair from the given reaction
	ClusterPair pair((PSICluster *) reaction->first,
			(PSICluster *) reaction->second);
	// Add the pair to the list
	reactingPairs.push_back(pair);
	// Setup the connectivity array
	setReactionConnectivity(reaction->first->getId());
	setReactionConnectivity(reaction->second->getId());

	return;
}

void PSICluster::createCombination(
		std::shared_ptr<ProductionReaction> reaction) {
	setReactionConnectivity(id);
	// Look for the other cluster
	IReactant * secondCluster;
	if (reaction->first->getId() == id)
		secondCluster = reaction->second;
	else
		secondCluster = reaction->first;

	// Creates the combining cluster
	CombiningCluster combCluster((PSICluster *) secondCluster);
	// Push the product into the list of clusters that combine with this one
	combiningReactants.push_back(combCluster);

	// Setup the connectivity array
	setReactionConnectivity(id);
	setReactionConnectivity(secondCluster->getId());

	return;
}

void PSICluster::createDissociation(
		std::shared_ptr<DissociationReaction> reaction) {
	// Look for the other cluster
	IReactant * emittedCluster;
	if (reaction->first->getId() == id)
		emittedCluster = reaction->second;
	else
		emittedCluster = reaction->first;

	// Create the pair of them where it is important that the
	// dissociating cluster is the first one
	ClusterPair pair((PSICluster *) reaction->dissociating,
			(PSICluster *) emittedCluster);
	// Add the pair to the dissociating pair vector
	dissociatingPairs.push_back(pair);

	// Setup the connectivity array
	setDissociationConnectivity(reaction->dissociating->getId());

	return;
}

void PSICluster::createEmission(
		std::shared_ptr<DissociationReaction> reaction) {
	// Create the pair of emitted clusters
	ClusterPair pair((PSICluster *) reaction->first,
			(PSICluster *) reaction->second);
	// Add the pair to the emission pair vector
	emissionPairs.push_back(pair);

	// Setup the connectivity array to itself
	setReactionConnectivity(id);

	return;
}

void PSICluster::optimizeReactions() {
	// Loop on the pairs to add reactions to the network
	for (auto it = reactingPairs.begin(); it != reactingPairs.end(); it++) {
		// Create the corresponding production reaction
		auto newReaction = std::make_shared<ProductionReaction>((*it).first,
				(*it).second);
		// Add it to the network
		newReaction = network->addProductionReaction(newReaction);
		// Link it to the pair
		(*it).reaction = newReaction;
	}

	for (auto it = combiningReactants.begin(); it != combiningReactants.end();
			it++) {
		// Create the corresponding production reaction
		auto newReaction = std::make_shared<ProductionReaction>((*it).combining,
				this);
		// Add it to the network
		newReaction = network->addProductionReaction(newReaction);
		// Link it to the pair
		(*it).reaction = newReaction;
	}

	for (auto it = dissociatingPairs.begin(); it != dissociatingPairs.end();
			it++) {
		// Create the corresponding dissociation reaction
		auto newReaction = std::make_shared<DissociationReaction>((*it).first,
				(*it).second, this);
		// Add it to the network
		newReaction = network->addDissociationReaction(newReaction);
		// Link it to the pair
		(*it).reaction = newReaction;
	}

	for (auto it = emissionPairs.begin(); it != emissionPairs.end(); it++) {
		// Create the corresponding dissociation reaction
		auto newReaction = std::make_shared<DissociationReaction>(this,
				(*it).first, (*it).second);
		// Add it to the network
		newReaction = network->addDissociationReaction(newReaction);
		// Link it to the pair
		(*it).reaction = newReaction;
	}

	return;
}

static std::vector<int> getFullConnectivityVector(std::set<int> connectivitySet,
		int size) {
	// Create a vector of zeroes with size equal to the network size
	std::vector<int> connectivity(size);

	// Set the value of the connectivity array to one for each element that is
	// in the set.
	for (auto it = connectivitySet.begin(); it != connectivitySet.end(); ++it) {
		connectivity[*it - 1] = 1;
	}

	return connectivity;
}

std::vector<int> PSICluster::getReactionConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(reactionConnectivitySet, network->getDOF());
}

std::vector<int> PSICluster::getDissociationConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(dissociationConnectivitySet,
			network->getDOF());
}

void PSICluster::resetConnectivities() {
	// Shrink the arrays to save some space. (About 10% or so.)
	reactingPairs.shrink_to_fit();
	combiningReactants.shrink_to_fit();
	dissociatingPairs.shrink_to_fit();
	emissionPairs.shrink_to_fit();

	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	setReactionConnectivity(heMomId);
	setDissociationConnectivity(heMomId);
	setReactionConnectivity(vMomId);
	setDissociationConnectivity(vMomId);

	// Loop on the effective reacting pairs
	for (auto it = reactingPairs.begin(); it != reactingPairs.end(); ++it) {
		// The cluster is connecting to both clusters in the pair
		setReactionConnectivity((*it).first->id);
		setReactionConnectivity((*it).second->id);
		setReactionConnectivity((*it).first->heMomId);
		setReactionConnectivity((*it).second->heMomId);
		setReactionConnectivity((*it).first->vMomId);
		setReactionConnectivity((*it).second->vMomId);
	}

	// Loop on the effective combining reactants
	for (auto it = combiningReactants.begin(); it != combiningReactants.end();
			++it) {
		// The cluster is connecting to the combining cluster
		setReactionConnectivity((*it).combining->id);
		setReactionConnectivity((*it).combining->heMomId);
		setReactionConnectivity((*it).combining->vMomId);
	}

	// Loop on the effective dissociating pairs
	for (auto it = dissociatingPairs.begin(); it != dissociatingPairs.end();
			++it) {
		// The cluster is connecting to the dissociating cluster which
		// is the first one by definition
		setDissociationConnectivity((*it).first->id);
		setDissociationConnectivity((*it).first->heMomId);
		setDissociationConnectivity((*it).first->vMomId);
	}

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	return;
}

void PSICluster::setReactionNetwork(
		const std::shared_ptr<IReactionNetwork> reactionNetwork) {
	// Call the superclass's method to actually set the reference
	Reactant::setReactionNetwork(reactionNetwork);

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	return;
}


double PSICluster::getDissociationFlux() const {
	// Initial declarations
	int nPairs = 0;
	double flux = 0.0;
	PSICluster *dissociatingCluster = nullptr;

	// Set the total number of reactants that dissociate to form this one
	nPairs = dissociatingPairs.size();
	// Loop over all dissociating clusters that form this cluster
	for (int j = 0; j < nPairs; j++) {
		// Get the dissociating cluster
		dissociatingCluster = dissociatingPairs[j].first;
		// Calculate the Dissociation flux
		flux += dissociatingPairs[j].reaction->kConstant
				* (double) dissociatingPairs[j].multiplicity
				* dissociatingCluster->getConcentration(
						dissociatingPairs[j].firstHeDistance,
						dissociatingPairs[j].firstVDistance);
	}

	// Return the flux
	return flux;
}

double PSICluster::getEmissionFlux() const {
	// Initial declarations
	int nPairs = 0;
	double flux = 0.0;

	// Set the total number of emission pairs
	nPairs = emissionPairs.size();
	// Loop over all the pairs
	for (int i = 0; i < nPairs; i++) {
		// Update the flux
		flux += emissionPairs[i].reaction->kConstant
				* (double) emissionPairs[i].multiplicity;
	}

	return flux * concentration;
}

double PSICluster::getProductionFlux() const {
	// Local declarations
	double flux = 0.0;
	int nPairs = 0;
	PSICluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Set the total number of reacting pairs
	nPairs = reactingPairs.size();
	// Loop over all the reacting pairs
	for (int i = 0; i < nPairs; i++) {
		// Get the two reacting clusters
		firstReactant = reactingPairs[i].first;
		secondReactant = reactingPairs[i].second;
		// Update the flux
		flux += reactingPairs[i].reaction->kConstant
				* (double) reactingPairs[i].multiplicity
				* firstReactant->getConcentration(
						reactingPairs[i].firstHeDistance,
						reactingPairs[i].firstVDistance)
				* secondReactant->getConcentration(
						reactingPairs[i].secondHeDistance,
						reactingPairs[i].secondVDistance);
	}

	// Return the production flux
	return flux;
}

double PSICluster::getCombinationFlux() const {
	// Local declarations
	double flux = 0.0;
	int nReactants = 0;
	PSICluster *combiningCluster = nullptr;

	// Set the total number of reactants that combine to form this one
	nReactants = combiningReactants.size();
	// Loop over all possible clusters
	for (int j = 0; j < nReactants; j++) {
		// Get the cluster that combines with this one
		combiningCluster = combiningReactants[j].combining;
		// Calculate the combination flux
		flux += combiningReactants[j].reaction->kConstant
				* (double) combiningReactants[j].multiplicity
				* combiningCluster->getConcentration(
						combiningReactants[j].heDistance,
						combiningReactants[j].vDistance);
	}

	return flux * concentration;
}

std::vector<double> PSICluster::getPartialDerivatives() const {
	// Local Declarations
	std::vector<double> partials(network->getDOF(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return partials;
}

void PSICluster::getPartialDerivatives(std::vector<double> & partials) const {
	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return;
}

void PSICluster::getProductionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numReactants = 0, index = 0;
	double value = 0.0;

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A
	numReactants = reactingPairs.size();
	for (int i = 0; i < numReactants; i++) {
		// Compute the contribution from the first part of the reacting pair
		value = reactingPairs[i].reaction->kConstant
				* (double) reactingPairs[i].multiplicity
				* reactingPairs[i].second->getConcentration(
						reactingPairs[i].secondHeDistance,
						reactingPairs[i].secondVDistance);
		index = reactingPairs[i].first->id - 1;
		partials[index] += value;
		index = reactingPairs[i].first->heMomId - 1;
		partials[index] += value * reactingPairs[i].firstHeDistance;
		index = reactingPairs[i].first->vMomId - 1;
		partials[index] += value * reactingPairs[i].firstVDistance;
		// Compute the contribution from the second part of the reacting pair
		value = reactingPairs[i].reaction->kConstant
				* (double) reactingPairs[i].multiplicity
				* reactingPairs[i].first->getConcentration(
						reactingPairs[i].firstHeDistance,
						reactingPairs[i].firstVDistance);
		index = reactingPairs[i].second->id - 1;
		partials[index] += value;
		index = reactingPairs[i].second->heMomId - 1;
		partials[index] += value * reactingPairs[i].secondHeDistance;
		index = reactingPairs[i].second->vMomId - 1;
		partials[index] += value * reactingPairs[i].secondVDistance;
	}

	return;
}

void PSICluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numReactants = 0, otherIndex = 0;
	PSICluster *cluster = nullptr;
	double value = 0.0;

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A
	numReactants = combiningReactants.size();
	for (int i = 0; i < numReactants; i++) {
		cluster = (PSICluster *) combiningReactants[i].combining;
		// Remember that the flux due to combinations is OUTGOING (-=)!
		// Compute the contribution from this cluster
		partials[id - 1] -= combiningReactants[i].reaction->kConstant
				* (double) combiningReactants[i].multiplicity
				* cluster->getConcentration(combiningReactants[i].heDistance,
						combiningReactants[i].vDistance);
		// Compute the contribution from the combining cluster
		value = combiningReactants[i].reaction->kConstant
				* (double) combiningReactants[i].multiplicity * concentration;
		otherIndex = cluster->id - 1;
		partials[otherIndex] -= value;
		otherIndex = cluster->heMomId - 1;
		partials[otherIndex] -= value * combiningReactants[i].heDistance;
		otherIndex = cluster->vMomId - 1;
		partials[otherIndex] -= value * combiningReactants[i].vDistance;
	}

	return;
}

void PSICluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numPairs = 0, index = 0;
	PSICluster *cluster = nullptr;
	double value = 0.0;

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)
	numPairs = dissociatingPairs.size();
	for (int i = 0; i < numPairs; i++) {
		// Get the dissociating cluster
		cluster = dissociatingPairs[i].first;
		value = dissociatingPairs[i].reaction->kConstant
				* (double) dissociatingPairs[i].multiplicity;
		index = cluster->id - 1;
		partials[index] += value;
		index = cluster->heMomId - 1;
		partials[index] += value * dissociatingPairs[i].firstHeDistance;
		index = cluster->vMomId - 1;
		partials[index] += value * dissociatingPairs[i].firstVDistance;
	}

	return;
}

void PSICluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numPairs = 0, index = 0;

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)
	numPairs = emissionPairs.size();
	for (int i = 0; i < numPairs; i++) {
		// Modify the partial derivative. Remember that the flux
		// due to emission is OUTGOING (-=)!
		index = id - 1;
		partials[index] -= emissionPairs[i].reaction->kConstant
				* (double) emissionPairs[i].multiplicity;
	}

	return;
}

void PSICluster::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

void PSICluster::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

double PSICluster::getLeftSideRate() const {
	// Initialize the rate and the cluster pointer
	double totalRate = 0.0;
	PSICluster *cluster = nullptr;

	// Loop on the combining reactants
	for (int i = 0; i < combiningReactants.size(); i++) {
		cluster = (PSICluster *) combiningReactants[i].combining;
		// Add the rate to the total rate
		totalRate += combiningReactants[i].reaction->kConstant
				* (double) combiningReactants[i].multiplicity
				* cluster->concentration;
	}

	// Loop on the emission pairs
	for (int i = 0; i < emissionPairs.size(); i++) {
		// Add the rate to the total rate
		totalRate += emissionPairs[i].reaction->kConstant
				* (double) emissionPairs[i].multiplicity;
	}

	return totalRate;
}

std::vector<int> PSICluster::getConnectivity() const {
	int connectivityLength = network->getDOF();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);
	auto reactionConnVector = getReactionConnectivity();
	auto dissociationConnVector = getDissociationConnectivity();

	// The reaction and dissociation vectors must have a length equal to the
	// number of clusters
	if (reactionConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The reaction vector is an incorrect length");
	}
	if (dissociationConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The dissociation vector is an incorrect length");
	}

	// Merge the two vectors such that the final vector contains
	// a 1 at a position if either of the connectivity arrays
	// have a 1
	for (int i = 0; i < connectivityLength; i++) {
		// Consider each connectivity array only if its type is enabled
		connectivity[i] = reactionConnVector[i] || dissociationConnVector[i];
	}

	return connectivity;
}
