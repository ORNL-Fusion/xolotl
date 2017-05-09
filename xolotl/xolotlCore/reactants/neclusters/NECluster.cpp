#include "NECluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

NECluster::NECluster() :
		Reactant() {
	// Set the reactant name appropriately
	name = "NECluster";

	// Setup the composition map.
	compositionMap[xeType] = 0;

	return;
}

NECluster::NECluster(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		Reactant(registry) {
	// Set the reactant name appropriately
	name = "NECluster";

	// Setup the composition map.
	compositionMap[xeType] = 0;

	return;
}

// The copy constructor
NECluster::NECluster(NECluster &other) :
		Reactant(other), reactingPairs(other.reactingPairs), combiningReactants(
				other.combiningReactants), dissociatingPairs(
				other.dissociatingPairs), emissionPairs(other.emissionPairs) {
	compositionMap[xeType] = other.compositionMap[xeType];

	// Recompute all of the temperature-dependent quantities
	setTemperature(other.getTemperature());

	return;
}

void NECluster::createProduction(std::shared_ptr<ProductionReaction> reaction) {
	// Create a cluster pair from the given reaction
	ClusterPair pair((NECluster *) reaction->first,
			(NECluster *) reaction->second);
	// Add the pair to the list
	reactingPairs.push_back(pair);
	// Setup the connectivity array
	setReactionConnectivity(reaction->first->getId());
	setReactionConnectivity(reaction->second->getId());

	return;
}

void NECluster::createCombination(
		std::shared_ptr<ProductionReaction> reaction) {
	setReactionConnectivity(id);
	// Look for the other cluster
	IReactant * secondCluster;
	if (reaction->first->getId() == id)
		secondCluster = reaction->second;
	else
		secondCluster = reaction->first;

	// Creates the combining cluster
	CombiningCluster combCluster((NECluster *) secondCluster);
	// Push the product into the list of clusters that combine with this one
	combiningReactants.push_back(combCluster);

	// Setup the connectivity array
	setReactionConnectivity(id);
	setReactionConnectivity(secondCluster->getId());

	return;
}

void NECluster::createDissociation(
		std::shared_ptr<DissociationReaction> reaction) {
	// Look for the other cluster
	IReactant * emittedCluster;
	if (reaction->first->getId() == id)
		emittedCluster = reaction->second;
	else
		emittedCluster = reaction->first;

	// Create the pair of them where it is important that the
	// dissociating cluster is the first one
	ClusterPair pair((NECluster *) reaction->dissociating,
			(NECluster *) emittedCluster);
	// Add the pair to the dissociating pair vector
	dissociatingPairs.push_back(pair);

	// Setup the connectivity array
	setDissociationConnectivity(reaction->dissociating->getId());

	return;
}

void NECluster::createEmission(std::shared_ptr<DissociationReaction> reaction) {
	// Create the pair of emitted clusters
	ClusterPair pair((NECluster *) reaction->first,
			(NECluster *) reaction->second);
	// Add the pair to the emission pair vector
	emissionPairs.push_back(pair);

	// Setup the connectivity array to itself
	setReactionConnectivity(id);

	return;
}

void NECluster::optimizeReactions() {
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
	for (auto it = combiningReactants.begin(); it != combiningReactants.end(); it++) {
		// Create the corresponding production reaction
		auto newReaction = std::make_shared<ProductionReaction>((*it).combining,
				this);
		// Add it to the network
		newReaction = network->addProductionReaction(newReaction);
		// Link it to the pair
		(*it).reaction = newReaction;
	}
	for (auto it = dissociatingPairs.begin(); it != dissociatingPairs.end(); it++) {
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
		auto newReaction = std::make_shared<DissociationReaction>(this, (*it).first,
				(*it).second);
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
		connectivity[(*it) - 1] = 1;
	}

	return connectivity;
}

std::vector<int> NECluster::getReactionConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(reactionConnectivitySet, network->getDOF());
}

std::vector<int> NECluster::getDissociationConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(dissociationConnectivitySet,
			network->getDOF());
}

void NECluster::resetConnectivities() {
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
	setReactionConnectivity(momId);
	setDissociationConnectivity(momId);

	// Loop on the reacting pairs
	for (auto it = reactingPairs.begin(); it != reactingPairs.end(); ++it) {
		// The cluster is connecting to both clusters in the pair
		setReactionConnectivity((*it).first->id);
		setReactionConnectivity((*it).first->momId);
		setReactionConnectivity((*it).second->id);
		setReactionConnectivity((*it).second->momId);
	}

	// Loop on the combining reactants
	for (auto it = combiningReactants.begin(); it != combiningReactants.end();
			++it) {
		// The cluster is connecting to the combining cluster
		setReactionConnectivity((*it).combining->id);
		setReactionConnectivity((*it).combining->momId);
	}

	// Loop on the effective dissociating pairs
	for (auto it = dissociatingPairs.begin(); it != dissociatingPairs.end();
			++it) {
		// The cluster is connecting to the dissociating cluster which
		// is the first one by definition
		setDissociationConnectivity((*it).first->id);
		setDissociationConnectivity((*it).first->momId);
	}

	// Don't loop on the emission pairs because
	// this cluster is not connected to them

	return;
}

void NECluster::setReactionNetwork(
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

double NECluster::getMoment() const {
	return 0.0;
}

double NECluster::getTotalFlux() {
	// Get the fluxes
	double prodFlux = getProductionFlux();
	double dissFlux = getDissociationFlux();
	double combFlux = getCombinationFlux();
	double emissFlux = getEmissionFlux();

	return prodFlux - combFlux + dissFlux - emissFlux;
}

double NECluster::getDissociationFlux() const {
	// Initial declarations
	int nPairs = 0;
	double flux = 0.0;
	NECluster *dissociatingCluster = nullptr;

	// Set the total number of reactants that dissociate to form this one
	nPairs = dissociatingPairs.size();
	// Loop over all dissociating clusters that form this cluster
	for (int j = 0; j < nPairs; j++) {
		// Get the dissociating cluster
		dissociatingCluster = dissociatingPairs[j].first;
		// Calculate the Dissociation flux
		flux += dissociatingPairs[j].reaction->kConstant
				* dissociatingCluster->getConcentration(
						dissociatingPairs[j].firstDistance);
	}

	// Return the flux
	return flux;
}

double NECluster::getEmissionFlux() const {
	// Initial declarations
	int nPairs = 0;
	double flux = 0.0;

	// Set the total number of emission pairs
	nPairs = emissionPairs.size();
	// Loop over all the pairs
	for (int i = 0; i < nPairs; i++) {
		// Update the flux
		flux += emissionPairs[i].reaction->kConstant;
	}

	return flux * concentration;
}

double NECluster::getProductionFlux() const {
	// Local declarations
	double flux = 0.0;
	int nPairs = 0;
	NECluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Set the total number of reacting pairs
	nPairs = reactingPairs.size();
	// Loop over all the reacting pairs
	for (int i = 0; i < nPairs; i++) {
		// Get the two reacting clusters
		firstReactant = reactingPairs[i].first;
		secondReactant = reactingPairs[i].second;
		// Update the flux
		flux += reactingPairs[i].reaction->kConstant
				* firstReactant->getConcentration(
						reactingPairs[i].firstDistance)
				* secondReactant->getConcentration(
						reactingPairs[i].secondDistance);
	}

	// Return the production flux
	return flux;
}

double NECluster::getCombinationFlux() const {
	// Local declarations
	double flux = 0.0;
	int nReactants = 0;
	NECluster *combiningCluster = nullptr;

	// Set the total number of reactants that combine to form this one
	nReactants = combiningReactants.size();
	// Loop over all possible clusters
	for (int j = 0; j < nReactants; j++) {
		// Get the cluster that combines with this one
		combiningCluster = combiningReactants[j].combining;
		// Calculate Second term of production flux
		flux += combiningReactants[j].reaction->kConstant
				* combiningCluster->getConcentration(
						combiningReactants[j].distance);
	}

	return flux * concentration;
}

std::vector<double> NECluster::getPartialDerivatives() const {
	// Local Declarations
	std::vector<double> partials(network->getDOF(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return partials;
}

void NECluster::getPartialDerivatives(std::vector<double> & partials) const {
	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return;
}

void NECluster::getProductionPartialDerivatives(
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
				* reactingPairs[i].second->getConcentration(
						reactingPairs[i].secondDistance);
		index = reactingPairs[i].first->id - 1;
		partials[index] += value;
		index = reactingPairs[i].first->momId - 1;
		partials[index] += value * reactingPairs[i].firstDistance;
		// Compute the contribution from the second part of the reacting pair
		value = reactingPairs[i].reaction->kConstant
				* reactingPairs[i].first->getConcentration(
						reactingPairs[i].firstDistance);
		index = reactingPairs[i].second->id - 1;
		partials[index] += value;
		index = reactingPairs[i].second->momId - 1;
		partials[index] += value * reactingPairs[i].secondDistance;
	}

	return;
}

void NECluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numReactants = 0, otherIndex = 0;
	NECluster *cluster = nullptr;
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
		cluster = (NECluster *) combiningReactants[i].combining;
		// Remember that the flux due to combinations is OUTGOING (-=)!
		// Compute the contribution from this cluster
		partials[id - 1] -= combiningReactants[i].reaction->kConstant
				* cluster->getConcentration(combiningReactants[i].distance);
		// Compute the contribution from the combining cluster
		value = combiningReactants[i].reaction->kConstant * concentration;
		// Get the index of cluster
		otherIndex = cluster->id - 1;
		partials[otherIndex] -= value;
		otherIndex = cluster->momId - 1;
		partials[otherIndex] -= value * combiningReactants[i].distance;
	}

	return;
}

void NECluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numPairs = 0, index = 0;
	NECluster *cluster = nullptr;

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
		index = cluster->id - 1;
		partials[index] += dissociatingPairs[i].reaction->kConstant;
		index = cluster->momId - 1;
		partials[index] += dissociatingPairs[i].reaction->kConstant
				* dissociatingPairs[i].firstDistance;
	}

	return;
}

void NECluster::getEmissionPartialDerivatives(
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
		partials[index] -= emissionPairs[i].reaction->kConstant;
	}

	return;
}

void NECluster::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

void NECluster::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

double NECluster::getLeftSideRate() const {
	// Initialize the rate and the cluster pointer
	double totalRate = 0.0;
	NECluster *cluster = nullptr;

	// Loop on the combining reactants
	for (int i = 0; i < combiningReactants.size(); i++) {
		cluster = (NECluster *) combiningReactants[i].combining;
		// Add the rate to the total rate
		totalRate += combiningReactants[i].reaction->kConstant
				* cluster->concentration;
	}

	// Loop on the emission pairs
	for (int i = 0; i < emissionPairs.size(); i++) {
		// Add the rate to the total rate
		totalRate += emissionPairs[i].reaction->kConstant;
	}

	return totalRate;
}

std::vector<int> NECluster::getConnectivity() const {
	int connectivityLength = network->getDOF();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);
	auto reactionConnVector = getReactionConnectivity();
	auto dissociationConnVector = getDissociationConnectivity();

	// The reaction and dissociation vectors must have a length equal to the
	// number of clusters
	if (reactionConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The reaction vector has an incorrect length");
	}
	if (dissociationConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The dissociation vector has an incorrect length");
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
