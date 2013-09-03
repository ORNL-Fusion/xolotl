#include "PSICluster.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

PSICluster::PSICluster() {
	// Set the size
	size = 1;
	// Zero out the binding energies
	bindingEnergies.resize(4);
	bindingEnergies[0] = 0.0;
	bindingEnergies[1] = 0.0;
	bindingEnergies[2] = 0.0;
	bindingEnergies[3] = 0.0;
	// Zero out the diffusion factor and migration energy
	diffusionFactor = 0.0;
	migrationEnergy = 0.0;
	// Set the reactant name appropriately
	name = "PSICluster";
}

PSICluster::PSICluster(const int clusterSize) :
		Reactant() {

	// Set the size
	size = (clusterSize > 0) ? clusterSize : 1;
	// Zero out the binding energies
	bindingEnergies.resize(4);
	bindingEnergies[0] = 0.0;
	bindingEnergies[1] = 0.0;
	bindingEnergies[2] = 0.0;
	bindingEnergies[3] = 0.0;
	// Zero out the diffusion factor and migration energy
	diffusionFactor = 0.0;
	migrationEnergy = 0.0;
	// Set the reactant name appropriately
	name = "PSICluster";
}

// The copy constructor with a huge initialization list!
PSICluster::PSICluster(const PSICluster &other) :
		Reactant(other), size(other.size), diffusionFactor(
				other.diffusionFactor), bindingEnergies(other.bindingEnergies), migrationEnergy(
				other.migrationEnergy), reactionConnectivity(
				other.reactionConnectivity), dissociationConnectivity(
				other.dissociationConnectivity), reactingPairs(
				other.reactingPairs), combiningReactants(
				other.combiningReactants) {
}

std::shared_ptr<Reactant> PSICluster::clone() {
	std::shared_ptr<Reactant> reactant(new PSICluster(*this));
	return reactant;
}

PSICluster::~PSICluster() {
}

int PSICluster::getSize() const {
	// Return this cluster's size
	return size;
}

void PSICluster::setReactionNetwork(
		const std::shared_ptr<ReactionNetwork> reactionNetwork) {

	// Call the superclass's method to actually set the reference
	Reactant::setReactionNetwork(reactionNetwork);

	// ----- Handle the connectivty for PSIClusters -----

	// Extract properties from the network
	auto properties = network->getProperties();
	int connectivityLength = network->getAll()->size();

	// Get the enabled reaction type flags
	bool reactionsEnabled = (properties["reactionsEnabled"] == "true");
	bool dissociationsEnabled = (properties["dissociationsEnabled"] == "true");

	// Resize the connectivity arrays
	reactionConnectivity.resize(connectivityLength, 0);
	dissociationConnectivity.resize(connectivityLength, 0);
	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();

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
}

double PSICluster::getDissociationFlux(double temperature) const {

	int nReactants = 0, oneIndex = -1;
	double diss = 0.0, conc = 0.0;
	std::shared_ptr<PSICluster> first, second;
	std::map<std::string, int> oneHe, oneV, oneI;

	// Only try this if the network is available
	if (network != NULL) {
		// Get the total number of reactants in the network
		nReactants = network->size();
		// Get this PSICluster or subclasses' cluster map
		auto composition = getComposition();
		// Get the number of species to determine if this
		// cluster is mixed or single
		int numSpecies = (composition["He"] > 0) + (composition["V"] > 0)
				+ (composition["I"] > 0);
		// If no species, throw error
		if (numSpecies == 0) {
			// Bad if we have no species
			throw std::string("Cluster map contains no species");
		} else if (numSpecies == 1) {
			// We know we are a single species,
			// but we need to know which one so we can
			// get the correct species He, V, or I to calculate
			// the dissociation constant.
			if (composition["He"]) {
				second = std::dynamic_pointer_cast<PSICluster>(network->get("He", 1));
			} else if (composition["V"]) {
				second = std::dynamic_pointer_cast<PSICluster>(network->get("V", 1));
			} else if (composition["I"]) {
				second = std::dynamic_pointer_cast<PSICluster>(network->get("I", 1));
			}
			// Loop over all reactants and see if we
			// have a dissociation connection
			auto reactants = network->getAll(); //Assumed it is ordered exactly!
			for (int i = 0; i < nReactants; i++) {
				// Only calculate if we are connected
				if (dissociationConnectivity.at(i) == 1) {
					first = std::dynamic_pointer_cast<PSICluster>(reactants->at(i));
					// Calculate the dissociation flux
					diss += calculateDissociationConstant(*first, *second,
							temperature)
							* first->getConcentration();
				}
			}
		} else if (numSpecies == 2) {
			throw std::string(
					"Mixed Species dissociation flux must be implemented by subclass.");
		}

	}

	// Return the flux
	return diss;
}

double PSICluster::getProductionFlux(double temperature) const {

	// Local declarations
	double flux = 0.0, kPlus = 0.0;
	int thisClusterIndex = 0;
	std::shared_ptr<PSICluster> firstReactant, secondReactant;
	int nPairs = 0;

	// Only try this if the network is available
	if (network != NULL) {
		// Set the total number of reacting pairs
		nPairs = reactingPairs.size();
		// Loop over all the reacting pairs
		for (int i = 0; i < nPairs; i++) {
			// Get the reactants
			firstReactant = reactingPairs.at(i).first;
			secondReactant = reactingPairs.at(i).second;
			// Update the flux
			flux += calculateReactionRateConstant(*firstReactant, *secondReactant,
					temperature) * firstReactant->getConcentration()
					* secondReactant->getConcentration();
		}

	}

	// Return the production flux
	return flux;
}

/**
 * This operation returns the total change in this cluster due to
 * the combination of this cluster with others.
 * @param temperature The temperature at which to calculate the flux
 * @return The flux due to this cluster combining with other clusters.
 */
double PSICluster::getCombinationFlux(double temperature) const {
	// Local declarations
	double flux = 0.0, kPlus = 0.0;
	int thisClusterIndex = 0;
	std::shared_ptr<PSICluster> outerReactant, thisReactant;
	int nReactants = 0;

	// Only try this if the network is available
	if (network != NULL) {
		// Set the total network nReactants
		nReactants = combiningReactants.size();
		// Loop over all possible clusters
		for (int j = 0; j < nReactants; j++) {
			outerReactant = std::dynamic_pointer_cast<PSICluster>(combiningReactants.at(j));
			// Calculate Second term of production flux
			flux += calculateReactionRateConstant(*this, *outerReactant,
					temperature) * outerReactant->getConcentration();
		}

	}

	// Return the production flux
	return (flux * getConcentration());
}

double PSICluster::getTotalFlux(double temperature) const {
	return getProductionFlux(temperature) - getCombinationFlux(temperature)
			+ getDissociationFlux(temperature);
}

double PSICluster::getDiffusionFactor() const {
	// Return the diffusion factor
	return diffusionFactor;
}

void PSICluster::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	return;
}

double PSICluster::getDiffusionCoefficient(double temperature) const{
	// Use the Arrhenius equation to compute the diffusion coefficient
	double k_b = xolotlCore::kBoltzmann;
	double kernel = -migrationEnergy / (k_b * temperature);
	return diffusionFactor * exp(kernel);
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
	return;
}

double PSICluster::calculateReactionRateConstant (
		const PSICluster & firstReactant,
		const PSICluster & secondReactant,
		const double temperature) const {

	// Get the reaction radii
	double r_first = firstReactant.getReactionRadius();
	double r_second = secondReactant.getReactionRadius();

	// Get the diffusion coefficients
	double firstDiffusion =  firstReactant.getDiffusionCoefficient(temperature);
	double secondDiffusion = secondReactant.getDiffusionCoefficient(temperature);

	// Calculate and return
	double k_plus = 4.0 * xolotlCore::pi * (r_first + r_second)
			* (firstDiffusion + secondDiffusion);
	return k_plus;
}

double PSICluster::calculateDissociationConstant(
		const PSICluster & firstReactant,
		const PSICluster & secondReactant,
		const double temperature) const {

	// Local Declarations
	int bindingEnergyIndex = -1;
	double ra = 1, rb = 1;
	double atomicVolume = 1.0;
	auto secondComposition = secondReactant.getComposition();

	// Get the binding energy index
	if (secondComposition["He"] == 1 && secondComposition["V"] == 0 && secondComposition["I"] == 0) {
		bindingEnergyIndex = 0;
	} else if (secondComposition["He"] == 0 && secondComposition["V"] == 1
			&& secondComposition["I"] == 0) {
		bindingEnergyIndex = 1;
	} else if (secondComposition["He"] == 0 && secondComposition["V"] == 0
			&& secondComposition["I"] == 1) {
		bindingEnergyIndex = 2;
	} else {
		return 0.0;
	}

	// Calculate the Reaction Rate Constant
	double kPlus = calculateReactionRateConstant(firstReactant, secondReactant,
			temperature);

	// Calculate and return
	double k_minus_exp = exp(
			bindingEnergies.at(bindingEnergyIndex)
					/ (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;
	return k_minus;
}

bool PSICluster::isProductReactant(const Reactant & reactantI,
		const Reactant & reactantJ) {
	// Base class should just return false
	return false;
}

const double PSICluster::getReactionRadius() const {
	return 0.0;
}

std::shared_ptr<std::vector<int>> PSICluster::getConnectivity() const {

	int connectivityLength = network->size();
	std::shared_ptr<std::vector<int>> connectivity = std::make_shared<std::vector<int>>(connectivityLength,0);

	// The reaction and dissociate vectors must be the same length
	// as the number of reactants
	if (reactionConnectivity.size() != connectivityLength) {
		throw std::string("The reaction vector is an incorrect length");
	}

	if (dissociationConnectivity.size() != connectivityLength) {
		throw std::string("The dissociation vector is an incorrect length");
	}

	// Merge the two vectors such that the final vector contains
	// a 1 at a position if either of the connectivity arrays
	// have a 1
	for (int i = 0; i < connectivityLength; i++) {
		// Consider each connectivity array only if its type is enabled
		(*connectivity)[i] = reactionConnectivity[i]
				|| dissociationConnectivity[i];
	}

	return connectivity;
}

void PSICluster::createReactionConnectivity() {
	// By default, generate an array with a zero for each reactant
	// in the network
	reactionConnectivity.clear();
	reactionConnectivity.resize(network->size(), 0);
}

void PSICluster::createDissociationConnectivity() {
	// By default, generate an array with a zero for each reactant
	// in the network
	dissociationConnectivity.clear();
	dissociationConnectivity.resize(network->size(), 0);
}

std::map<std::string, int> PSICluster::getClusterMap() {
	// Create an empty cluster map
	std::map<std::string, int> clusterMap;
	return clusterMap;
}
