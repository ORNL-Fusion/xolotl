#include "PSICluster.h"
#include "HandlerRegistryFactory.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

// Create the static map of binding energies
// It has to be initialized here because it was not working otherwise
std::unordered_map<std::string, int> PSICluster::bindingEnergyIndexMap
	= {{heType, 0}, {vType, 1}, {iType, 2}};

PSICluster::PSICluster() :
		Reactant() {
	// Set the size
	size = 1;
	// Zero out the binding energies
	bindingEnergies.resize(3);
	bindingEnergies[0] = 0.0;
	bindingEnergies[1] = 0.0;
	bindingEnergies[2] = 0.0;
	// Zero out the diffusion factor and migration energy
	diffusionFactor = 0.0;
	diffusionCoefficient = 0.0;
	migrationEnergy = 0.0;
	// Set the reactant name appropriately
	name = "PSICluster";
	// No type name to set. Leave it empty.
	// Setup the composition map.
	compositionMap[heType] = 0;
	compositionMap[vType] = 0;
	compositionMap[iType] = 0;
	// Set the default reaction radius to 0. (Doesn't react.)
	reactionRadius = 0.0;
}

PSICluster::PSICluster(const int clusterSize,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		Reactant(registry) {

	// Set the size
	size = (clusterSize > 0) ? clusterSize : 1;
	// Zero out the binding energies
	bindingEnergies.resize(3);
	bindingEnergies[0] = 0.0;
	bindingEnergies[1] = 0.0;
	bindingEnergies[2] = 0.0;
	// Zero out the diffusion factor and migration energy
	diffusionFactor = 0.0;
	diffusionCoefficient = 0.0;
	migrationEnergy = 0.0;
	// Set the reactant name appropriately
	name = "PSICluster";
	// Setup the composition map.
	compositionMap[heType] = 0;
	compositionMap[vType] = 0;
	compositionMap[iType] = 0;
	// Set the default reaction radius to 0. (Doesn't react.)
	reactionRadius = 0.0;

	// Set up an event counter to count the number of times getDissociationFlux is called
	getDissociationFluxCounter = handlerRegistry->getEventCounter(
			"PSICluster_getDissociationFlux_Counter");
}

// The copy constructor with a huge initialization list!
PSICluster::PSICluster(const PSICluster &other) :
		Reactant(other), size(other.size), diffusionFactor(
				other.diffusionFactor), thisNetworkIndex(
				other.thisNetworkIndex), bindingEnergies(
				other.bindingEnergies), migrationEnergy(
				other.migrationEnergy), reactionRadius(
				other.reactionRadius), reactingPairs(
				other.reactingPairs), combiningReactants(
				other.combiningReactants), dissociatingPairs(
				other.dissociatingPairs), emissionPairs(
				other.emissionPairs), reactionConnectivitySet(
				other.reactionConnectivitySet), dissociationConnectivitySet(
				other.dissociationConnectivitySet) {

	// Recompute all of the temperature-dependent quantities
	setTemperature(other.getTemperature());

	// Set the single species cluster reference
	sameTypeSizeOneCluster = other.sameTypeSizeOneCluster;

	// Set up an event counter to count the number of times getDissociationFlux is called
	getDissociationFluxCounter = handlerRegistry->getEventCounter(
			"getDissociationFlux_Counter");

}

std::shared_ptr<Reactant> PSICluster::clone() {
	std::shared_ptr<Reactant> reactant(new PSICluster(*this));
	return reactant;
}

PSICluster::~PSICluster() {
}

void PSICluster::printReaction(const Reactant & firstReactant,
		const Reactant & secondReactant,
		const Reactant & productReactant) const {

	auto firstComp = firstReactant.getComposition();
	auto secondComp = secondReactant.getComposition();
	auto productComp = productReactant.getComposition();

	std::cout << firstReactant.getName() << "(" << firstComp[heType] << ", "
			<< firstComp[vType] << ", " << firstComp[iType] << ") + "
			<< secondReactant.getName() << "(" << secondComp[heType] << ", "
			<< secondComp[vType] << ", " << secondComp[iType] << ") -> "
			<< productReactant.getName() << "(" << productComp[heType] << ", "
			<< productComp[vType] << ", " << productComp[iType] << ")"
			<< std::endl;

	return;
}

void PSICluster::printDissociation(const Reactant & firstReactant,
		const Reactant & secondReactant,
		const Reactant & productReactant) const {

	auto firstComp = firstReactant.getComposition();
	auto secondComp = secondReactant.getComposition();
	auto productComp = productReactant.getComposition();

	std::cout << firstReactant.getName() << "(" << firstComp[heType] << ", "
			<< firstComp[vType] << ", " << firstComp[iType] << ") -> "
			<< secondReactant.getName() << "(" << secondComp[heType] << ", "
			<< secondComp[vType] << ", " << secondComp[iType] << ") + "
			<< productReactant.getName() << "(" << productComp[heType] << ", "
			<< productComp[vType] << ", " << productComp[iType] << ")"
			<< std::endl;

	return;
}

static std::vector<int> getFullConnectivityVector(std::set<int> connectivitySet,
		int size) {

	// Create a vector of zeroes with size equal to the network size
	std::vector<int> connectivity(size);

	// Set the value of the connectivity array to one for each element that is
	// in the set.
	for (auto it = connectivitySet.begin(); it != connectivitySet.end(); it++) {
		connectivity[*it - 1] = 1;
	}

	return connectivity;
}

void PSICluster::setReactionConnectivity(int clusterId) {
	// Add the cluster to the set.
	reactionConnectivitySet.insert(clusterId);
}

std::vector<int> PSICluster::getReactionConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(reactionConnectivitySet, network->size());
}

std::set<int> PSICluster::getReactionConnectivitySet() const {
	return std::set<int>(reactionConnectivitySet);
}

void PSICluster::setDissociationConnectivity(int clusterId) {
	// Add the cluster to the set.
	dissociationConnectivitySet.insert(clusterId);
}

std::vector<int> PSICluster::getDissociationConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(dissociationConnectivitySet,
			network->size());
}

const std::set<int> & PSICluster::getDissociationConnectivitySet() const {
	return dissociationConnectivitySet;
}

int PSICluster::getSize() const {
	// Return this cluster's size
	return size;
}

void PSICluster::setReactionNetwork(
		const std::shared_ptr<ReactionNetwork> reactionNetwork) {

	// Call the superclass's method to actually set the reference
	Reactant::setReactionNetwork(reactionNetwork);

	// Extract properties from the network
	auto properties = network->getProperties();
	int connectivityLength = network->size();

	// Get the enabled reaction type flags
	bool reactionsEnabled = (properties["reactionsEnabled"] == "true");
	bool dissociationsEnabled = (properties["dissociationsEnabled"] == "true");

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	// Get the index/id of this cluster in the reaction network.
	thisNetworkIndex = getId() - 1;

	// ----- Handle the connectivty for PSIClusters -----

	// Generate the reactant and dissociation connectivity arrays.
	// This only must be done once since the arrays are stored as
	// member attributes. Only perform these tasks if the reaction
	// types are enabled.
	if (reactionsEnabled) {
		createReactionConnectivity();
	}
	if (dissociationsEnabled) {
		createDissociationConnectivity();
	}

	// Shrink the arrays to save some space. (About 10% or so.)
	reactingPairs.shrink_to_fit();
	combiningReactants.shrink_to_fit();
	dissociatingPairs.shrink_to_fit();
	emissionPairs.shrink_to_fit();

	// Get the cluster that represents the same type as this cluster, but with
	// size = 1. This only works for single species clusters to being with.
	if (!isMixed()) {
		sameTypeSizeOneCluster = (PSICluster *) network->get(typeName, 1);
	}

	return;
}

double PSICluster::getDissociationFlux(double temperature) const {
	// increment the getDissociationFlux counter
	getDissociationFluxCounter->increment();

	// Initial declarations
	int nPairs = 0;
	double flux = 0.0, fluxMultiplier = 1.0;

	// Only try this if the network is available
	if (network != NULL) {
		// Set the total number of reactants that dissociate to form this one
		nPairs = dissociatingPairs.size();
		// Loop over all dissociating clusters that form this cluster
		for (int j = 0; j < nPairs; j++) {
			auto dissociatingCluster = dissociatingPairs[j].first;
			// The second element of the pair is the cluster that is also
			// emitted by the dissociation
			auto otherEmittedCluster = dissociatingPairs[j].second;
			// Calculate the Dissociation flux
			flux += fluxMultiplier
					* dissociatingPairs[j].kConstant
					* dissociatingCluster->getConcentration();

			// Need to be added twice when a cluster is emitted with itself
			// because it is just once in the dissociation pairs list
			if (this->getId() == otherEmittedCluster->getId()) {
				flux += fluxMultiplier
						* dissociatingPairs[j].kConstant
						* dissociatingCluster->getConcentration();
			}
		}
	}

	// Return the flux
	return flux;
}

double PSICluster::getEmissionFlux(double temperature) const {
	// Initial declarations
	int nPairs = 0;
	double flux = 0.0, fluxMultiplier = 1.0;

	// Only try this if the network is available
	if (network != NULL) {
		// Set the total number of emission pairs
		nPairs = emissionPairs.size();
		// Loop over all the pairs
		for (int i = 0; i < nPairs; i++) {
			auto firstCluster = emissionPairs[i].first;
			auto secondCluster = emissionPairs[i].second;
			// Update the flux
			flux += emissionPairs[i].kConstant;
		}
	}

	// Return the flux
	return (flux * getConcentration());
}

double PSICluster::getProductionFlux(double temperature) const {
	// Local declarations
	double flux = 0.0;
	double conc1 = 0.0, conc2 = 0.0;
	int nPairs = 0;

	// Only try this if the network is available
	if (network != NULL) {
		// Set the total number of reacting pairs
		nPairs = reactingPairs.size();
		// Loop over all the reacting pairs
		for (int i = 0; i < nPairs; i++) {
			// Get the reactants
			auto firstReactant = reactingPairs[i].first;
			auto secondReactant = reactingPairs[i].second;
			// Update the flux
			conc1 = firstReactant->getConcentration();
			conc2 = secondReactant->getConcentration();
			flux += reactingPairs[i].kConstant * conc1 * conc2;
		}
	}

	// Return the production flux
	return flux;
}

double PSICluster::getCombinationFlux(double temperature) const {

	// Local declarations
	double flux = 0.0, conc = 0.0;
	int nReactants = 0;

	// Set the total number of reactants that combine to form this one
	nReactants = combiningReactants.size();
	// Loop over all possible clusters
	for (int j = 0; j < nReactants; j++) {
		auto otherCluster = (PSICluster *) combiningReactants[j].combining;
		conc = otherCluster->getConcentration();
		// Calculate Second term of production flux
		flux += combiningReactants[j].kConstant * conc;

		// Need to be added twice when a cluster combine with itself
		// because it is just once in the combining reactant list
		if (this->getId() == otherCluster->getId()) {
			flux += combiningReactants[j].kConstant * conc;
		}
	}

	// Return the production flux
	return (flux * getConcentration());
}

double PSICluster::getTotalFlux(const double temperature) {

	// Get the fluxes
	double prodFlux, combFlux, dissFlux, emissFlux;
	prodFlux = getProductionFlux(temperature);
	combFlux = getCombinationFlux(temperature);
	dissFlux = getDissociationFlux(temperature);
	emissFlux = getEmissionFlux(temperature);

	return prodFlux - combFlux + dissFlux - emissFlux;
}

double PSICluster::getDiffusionFactor() const {
	// Return the diffusion factor
	return diffusionFactor;
}

void PSICluster::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(getTemperature());
	return;
}

double PSICluster::getDiffusionCoefficient() const {
	return diffusionCoefficient;
}

std::vector<double> PSICluster::getBindingEnergies() const {
	// Local Declarations
	std::vector<double> energyVector;

	// Set the return binding energies
	energyVector = bindingEnergies;

	// Return the energies
	return energyVector;
}

void PSICluster::setBindingEnergies(const std::vector<double> energies) {
	// Set the binding energies
	bindingEnergies = energies;
	return;
}

double PSICluster::getMigrationEnergy() const {
	// Return the migration energy
	return migrationEnergy;
}

void PSICluster::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(getTemperature());
	return;
}

double PSICluster::calculateReactionRateConstant(
		const PSICluster & firstReactant, const PSICluster & secondReactant,
		double temperature) const {

	// Get the reaction radii
	double r_first = firstReactant.getReactionRadius();
	double r_second = secondReactant.getReactionRadius();

	// Get the diffusion coefficients
	double firstDiffusion = firstReactant.getDiffusionCoefficient();
	double secondDiffusion = secondReactant.getDiffusionCoefficient();

	// Calculate and return
	double k_plus = 4.0 * xolotlCore::pi * (r_first + r_second)
			* (firstDiffusion + secondDiffusion);
	return k_plus;
}

double PSICluster::calculateDissociationConstant(const PSICluster & dissociatingCluster,
		const PSICluster & singleCluster, const PSICluster & secondCluster,
		double temperature) const {

	// Local Declarations
	int bindingEnergyIndex = -1;
	double atomicVolume = 1.0; // Currently calculated below for He and W only!
	double bindingEnergy = 0.0;

	// Get the binding energy index.
	bindingEnergyIndex = bindingEnergyIndexMap[singleCluster.typeName];
	// The atomic volume is computed by considering the BCC structure of the
	// tungsten. In a given lattice cell in tungsten there are tungsten atoms
	// at each corner and a tungsten atom in the center. The tungsten atoms at
	// the corners are shared across a total of eight cells. The fraction of
	// the volume of the lattice cell that is filled with tungsten atoms is the
	// atomic volume and is a_0^3/(8*1/8 + 1) = 0.5*a_0^3.
	atomicVolume = 0.5 * xolotlCore::latticeConstant
			* xolotlCore::latticeConstant * xolotlCore::latticeConstant;

	// Calculate the Reaction Rate Constant
	double kPlus = 0.0;
	kPlus = calculateReactionRateConstant(singleCluster, secondCluster,
			temperature);

	// Calculate and return
	bindingEnergy = dissociatingCluster.getBindingEnergies()[bindingEnergyIndex];
	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;
	return k_minus;
}

void PSICluster::computeRateConstants(double temperature) {
	// Compute the reaction constant associated to the reacting pairs
	// Set the total number of reacting pairs
	int nPairs = reactingPairs.size();
	// Loop on them
	for (int i = 0; i < nPairs; i++) {
		// Get the reactants
		auto firstReactant = reactingPairs[i].first;
		auto secondReactant = reactingPairs[i].second;
		// Compute the reaction constant
		auto rate = calculateReactionRateConstant(*firstReactant,
				*secondReactant, temperature);
		// Set it in the pair
		reactingPairs[i].kConstant = rate;
	}

	// Compute the reaction constant associated to the combining reactants
	// Set the total number of combining reactants
	int nReactants = combiningReactants.size();
	// Loop on them
	for (int i = 0; i < nReactants; i++) {
		// Get the reactants
		auto combiningReactant = combiningReactants[i].combining;
		// Compute the reaction constant
		auto rate = calculateReactionRateConstant(*this,
				*combiningReactant, temperature);
		// Set it in the combining cluster
		combiningReactants[i].kConstant = rate;
	}

	// Compute the dissociation constant associated to the dissociating clusters
	// Set the total number of dissociating clusters
	nPairs = dissociatingPairs.size();
	// Loop on them
	for (int i = 0; i < nPairs; i++) {
		auto dissociatingCluster = dissociatingPairs[i].first;
		// The second element of the pair is the cluster that is also
		// emitted by the dissociation
		auto otherEmittedCluster = dissociatingPairs[i].second;
		// Compute the dissociation constant
		double rate = 0.0;
		// The order of the cluster is important here because of the binding
		// energy used in the computation. It is taken from the type of the first cluster
		// which must be the single one
		if (size == 1) {
			// "this" is the single size one
			rate = calculateDissociationConstant(*dissociatingCluster,
					*this, *otherEmittedCluster, temperature);
		}
		else {
			// otherEmittedCluster is the single size one
			rate = calculateDissociationConstant(*dissociatingCluster,
					*otherEmittedCluster, *this, temperature);

		}
		// Set it in the pair
		dissociatingPairs[i].kConstant = rate;
	}

	// Compute the dissociation constant associated to the emission of pairs of clusters
	// Set the total number of emission pairs
	nPairs = emissionPairs.size();
	// Loop on them
	for (int i = 0; i < nPairs; i++) {
		auto firstCluster = emissionPairs[i].first;
		auto secondCluster = emissionPairs[i].second;
		// Compute the dissociation rate
		auto rate = calculateDissociationConstant(*this, *firstCluster,
				*secondCluster, temperature);
		// Set it in the pair
		emissionPairs[i].kConstant = rate;
	}

	return;
}

double PSICluster::getReactionRadius() const {
	return reactionRadius; // Computed by subclasses in constructors.
}

std::vector<int> PSICluster::getConnectivity() const {

	int connectivityLength = network->size();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);
	auto reactionConnVector = getReactionConnectivity();
	auto dissociationConnVector = getDissociationConnectivity();

	// The reaction and dissociation vectors must have a length equal to the
	// number of clusters
	if (reactionConnVector.size() != connectivityLength) {
		throw std::string("The reaction vector is an incorrect length");
	}
	if (dissociationConnVector.size() != connectivityLength) {
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

void PSICluster::createReactionConnectivity() {
	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(getId());

	// This cluster is always X_a

	// Initial declarations
	int firstSize = 0, secondSize = 0;

	// Single species clustering
	// X_(a-i) + X_i --> X_a
	for (firstSize = 1; firstSize <= (int) size / 2; firstSize++) {
		// Set the size of the second reactant
		secondSize = size - firstSize;
		// Get the first and second reactants for the reaction
		auto firstReactant = (PSICluster *) network->get(typeName, firstSize);
		auto secondReactant = (PSICluster *) network->get(typeName, secondSize);
		// Create a ReactingPair with the two reactants
		if (firstReactant && secondReactant) {
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

	// X_a + X_b --> X_(a+b)
	auto reactants = network->getAll(typeName);
	// combineClusters handles everything for this type of reaction
	combineClusters(reactants, typeName);

	return;
}

void PSICluster::createDissociationConnectivity() {

	// This cluster is always X_a

	// X_a --> X_(a-1) + X
	auto smallerReactant = network->get(typeName, size - 1);
	auto singleCluster = network->get(typeName, 1);
	emitClusters(singleCluster, smallerReactant);

	// X_(a+1) --> X_a + X
	auto biggerReactant = network->get(typeName, size + 1);
	dissociateCluster(biggerReactant, singleCluster);

	// Specific case for the single size cluster
	// for a = 1
	if (size == 1) {
		// all the cluster of the same type dissociate into it
		auto allSameTypeReactants = network->getAll(typeName);
		for (int i = 0; i < allSameTypeReactants.size(); i++) {
			// the one with size two was already added
			auto cluster = (PSICluster *) allSameTypeReactants[i];
			if (cluster->getSize() < 3)
				continue;

			// X_b is the dissociating one, X_(b-a) is the one
			// that is also emitted during the dissociation
			smallerReactant = network->get(typeName, cluster->getSize() - 1);
			dissociateCluster(cluster, smallerReactant);
		}
	}

	return;
}

void PSICluster::dissociateCluster(Reactant * dissociatingCluster,
		Reactant * emittedCluster) {
	// Test if the dissociatingCluster and the emittedCluster exist
	if (dissociatingCluster && emittedCluster) {
		// Cast to PSICluster so that we can get the information we need
		auto castedDissociatingCluster = (PSICluster *) dissociatingCluster;
		auto castedEmittedCluster = (PSICluster *) emittedCluster;

		// Create the pair of them where it is important that the
		// dissociating cluster is the first one
		// The dissociating constant will be computed later, it is set to 0.0 for now
		ClusterPair pair(castedDissociatingCluster, castedEmittedCluster, 0.0);
		// Add the pair to the dissociating pair vector
		// The connectivity is handled in emitCluster
		dissociatingPairs.push_back(pair);

		// Take care of the connectivity
		setDissociationConnectivity(dissociatingCluster->getId());
	}

	return;
}

void PSICluster::emitClusters(Reactant * firstEmittedCluster,
		Reactant * secondEmittedCluster) {
	// Test if the emitted clusters exist
	if (firstEmittedCluster && secondEmittedCluster) {
		// Cast to PSICluster so that we can get the information we need
		auto castedFirstCluster = (PSICluster *) firstEmittedCluster;
		auto castedSecondCluster = (PSICluster *) secondEmittedCluster;

		// Connect this cluster to itself since any reaction will affect it
		setDissociationConnectivity(getId());

		// Add the pair of emitted clusters to the vector of emissionPairs
		// The first cluster is the size one one
		// The dissociating constant will be computed later, it is set to 0.0 for now
		ClusterPair pair(castedFirstCluster, castedSecondCluster, 0.0);
		emissionPairs.push_back(pair);
	}

	return;
}

void PSICluster::getProductionPartialDerivatives(std::vector<double> & partials,
		double temperature) const {

	// Initial declarations
	int numReactants = 0, index = 0;
	double rateConstant = 0.0;

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
		index = reactingPairs[i].first->getId() - 1;
		partials[index] += reactingPairs[i].kConstant
				* reactingPairs[i].second->getConcentration();
		// Compute the contribution from the second part of the reacting pair
		index = reactingPairs[i].second->getId() - 1;
		partials[index] += reactingPairs[i].kConstant
				* reactingPairs[i].first->getConcentration();
	}

	return;
}

void PSICluster::getCombinationPartialDerivatives(
		std::vector<double> & partials, double temperature) const {

	// Initial declarations
	int numReactants = 0, otherIndex = 0;

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A
	numReactants = combiningReactants.size();
	for (int i = 0; i < numReactants; i++) {
		auto cluster = (PSICluster *) combiningReactants[i].combining;
		// Get the index of cluster
		otherIndex = cluster->getId() - 1;
		// Remember that the flux due to combinations is OUTGOING (-=)!
		// Compute the contribution from this cluster
		partials[thisNetworkIndex] -= combiningReactants[i].kConstant * cluster->getConcentration();
		// Compute the contribution from the combining cluster
		partials[otherIndex] -= combiningReactants[i].kConstant * getConcentration();

		// Need to be added twice when a cluster combine with itself
		// because it is just once in the combining reactant list
		if (this->getId() == cluster->getId()) {
			// Compute the contribution from this cluster
			partials[thisNetworkIndex] -= combiningReactants[i].kConstant * cluster->getConcentration();
			// Compute the contribution from the combining cluster
			partials[otherIndex] -= combiningReactants[i].kConstant * getConcentration();
		}
	}

	return;
}

void PSICluster::getDissociationPartialDerivatives(
		std::vector<double> & partials, double temperature) const {
	// Initial declarations
	int numPairs = 0, index = 0;

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)
	numPairs = dissociatingPairs.size();
	for (int i = 0; i < numPairs; i++) {
		// Get the dissociating cluster
		auto cluster = dissociatingPairs[i].first;
		auto emittedCluster = dissociatingPairs[i].second;
		index = cluster->getId() - 1;
		partials[index] += dissociatingPairs[i].kConstant;

		// Need to be added twice when a cluster is emitted with itself
		// because it is just once in the dissociation pairs list
		if (this->getId() == emittedCluster->getId()) {
			partials[index] += dissociatingPairs[i].kConstant;
		}
	}

	return;
}

void PSICluster::getEmissionPartialDerivatives(std::vector<double> & partials,
		double temperature) const {
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
		// Get the pair of clusters
		auto firstCluster = emissionPairs[i].first;
		auto secondCluster = emissionPairs[i].second;

		// Modify the partial derivative. Remember that the flux
		// due to emission is OUTGOING (-=)!
		index = getId() - 1;
		partials[index] -= emissionPairs[i].kConstant;

	}

	return;
}

std::vector<double> PSICluster::getPartialDerivatives(
		double temperature) const {
	// Local Declarations
	std::vector<double> partials(network->size(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials, temperature);
	getCombinationPartialDerivatives(partials, temperature);
	getDissociationPartialDerivatives(partials, temperature);
	getEmissionPartialDerivatives(partials, temperature);

	return partials;
}

void PSICluster::getPartialDerivatives(double temperature,
		std::vector<double> & partials) const {

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials, temperature);
	getCombinationPartialDerivatives(partials, temperature);
	getDissociationPartialDerivatives(partials, temperature);
	getEmissionPartialDerivatives(partials, temperature);

	return;
}

void PSICluster::combineClusters(std::vector<Reactant *> & reactants,
		std::string productName) {
	// Initial declarations
	std::map<std::string, int> myComposition = getComposition(),
			secondComposition;
	int numHe, numV, numI, secondNumHe, secondNumV, secondNumI;
	int otherId, productId;
	std::vector<int> compositionSizes { 0, 0, 0 };
	PSICluster * productCluster;
	// Setup the composition variables for this cluster
	numHe = myComposition[heType];
	numV = myComposition[vType];
	numI = myComposition[iType];

	int reactantVecSize = reactants.size();
	for (int i = 0; i < reactantVecSize; i++) {
		// Get the second reactant, its composition and its index
		auto secondCluster = (PSICluster *) reactants[i];
		secondComposition = secondCluster->getComposition();
		secondNumHe = secondComposition[heType];
		secondNumV = secondComposition[vType];
		secondNumI = secondComposition[iType];
		otherId = secondCluster->getId();
		int productSize = size + secondCluster->getSize();
		// Get and handle product for compounds
		if (productName == heVType || productName == heIType) {
			// Modify the composition vector
			compositionSizes[0] = numHe + secondNumHe;
			compositionSizes[1] = numV + secondNumV;
			compositionSizes[2] = numI + secondNumI;
			// Get the product
			productCluster = (PSICluster *) network->getCompound(productName,
					compositionSizes);
		} else {
			// Just get the product if it is a single-species
			productCluster = (PSICluster *) network->get(productName,
					productSize);
		}
		// React if the product exists in the network
		if (productCluster) {
			// Setup the connectivity array for the second reactant
			setReactionConnectivity(otherId);
			// Creates the combining cluster
			// The reaction constant will be computed later and is set to 0.0 for now
			CombiningCluster combCluster(secondCluster, 0.0);
			// Push the product onto the list of clusters that combine with this one
			combiningReactants.push_back(combCluster);
		}
	}

	return;
}

void PSICluster::replaceInCompound(std::vector<Reactant *> & reactants,
		std::string oldComponentName, std::string newComponentName) {
	// Local Declarations
	std::map<std::string, int> secondReactantComp,
			productReactantComp;
	int numReactants = reactants.size();
	int secondId = 0, productId = 0;

	// Loop over all of the extra reactants in this reaction and handle the replacement
	for (int i = 0; i < numReactants; i++) {
		// Get the second reactant and its composition
		auto secondReactant = (PSICluster *) reactants[i];
		secondReactantComp = secondReactant->getComposition();
		// Create the composition vector
		productReactantComp = secondReactantComp;
		// Updated the modified components
		productReactantComp[oldComponentName] =
				secondReactantComp[oldComponentName] - size;
		// Create the composition vector -- FIXME! This should be general!
		std::vector<int> productCompositionVector = { productReactantComp[heType],
				productReactantComp[vType], productReactantComp[iType] };
		// Get the product of the same type as the second reactant
		auto productReactant = network->getCompound(secondReactant->getType(),
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

void PSICluster::fillVWithI(std::string secondClusterName,
		std::vector<Reactant *> & reactants) {
	// Local Declarations
	std::shared_ptr<PSICluster> productCluster, thisCluster;
	std::string firstClusterName = typeName, productClusterName;
	int firstClusterSize = 0, secondClusterSize = 0, productClusterSize = 0;
	int secondId = 0, productId = 0, reactantVecSize = 0;

	// Get the number of V or I in this cluster (the "first")
	firstClusterSize = getSize();
	// Look at all of the second clusters, either V or I, and determine
	// if a connection exists.
	reactantVecSize = reactants.size();
	for (int i = 0; i < reactantVecSize; i++) {
		// Get the second cluster its size
		auto secondCluster = (PSICluster *) reactants[i];
		secondClusterSize = (secondCluster->getSize());
		// The only way this reaction is allowed is if the sizes are not equal.
		if (firstClusterSize != secondClusterSize) {
			// We have to switch on cluster type to make sure that the annihilation
			// is computed correctly.
			if (firstClusterName == vType) {
				// Compute the product size and set the product name for the case
				// where I is the second cluster
				if (secondClusterSize > firstClusterSize) {
					productClusterSize = secondClusterSize - firstClusterSize;
					productClusterName = iType;
				} else if (secondClusterSize < firstClusterSize) {
					productClusterSize = firstClusterSize - secondClusterSize;
					productClusterName = vType;
				}
			} else if (firstClusterName == iType) {
				// Compute the product size and set the product name for the case
				// where V is the second cluster
				if (firstClusterSize > secondClusterSize) {
					productClusterSize = firstClusterSize - secondClusterSize;
					productClusterName = iType;
				} else if (firstClusterSize < secondClusterSize) {
					productClusterSize = secondClusterSize - firstClusterSize;
					productClusterName = vType;
				}
			}
			// Get the product
			auto productCluster = (PSICluster *) network->get(
					productClusterName, productClusterSize);
			// Only deal with this reaction if the product exists. Otherwise the
			// whole reaction is forbidden.
			if (productCluster) {
				// Setup the connectivity array to handle the second reactant
				secondId = secondCluster->getId();
				setReactionConnectivity(secondId);
				// Creates the combining cluster
				// The reaction constant will be computed later and is set to 0.0 for now
				CombiningCluster combCluster(secondCluster, 0.0);
				// Push the second cluster onto the list of clusters that combine
				// with this one
				combiningReactants.push_back(combCluster);
			}
		}
	}
}

void PSICluster::recomputeDiffusionCoefficient(double temp) {

	// Return zero if the diffusion factor is zero.
	if (diffusionFactor == 0.0) {
		diffusionCoefficient = 0.0;
	} else {
		// Otherwise use the Arrhenius equation to compute the diffusion
		// coefficient
		double k_b = xolotlCore::kBoltzmann;
		double kernel = -1.0 * migrationEnergy / (k_b * temp);
		diffusionCoefficient = diffusionFactor * exp(kernel);
	}

	return;
}

void PSICluster::setTemperature(double temp) {
	// Set the temperature
	Reactant::setTemperature(temp);

	// Recompute the diffusion coefficient
	recomputeDiffusionCoefficient(temp);

	return;
}
