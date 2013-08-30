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
	// Create a smart pointer to a new connectivity array
	connectivity = std::shared_ptr<std::vector<int>>(new std::vector<int>);
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
	// Create a smart pointer to a new connectivity array
	connectivity = std::shared_ptr<std::vector<int>>(new std::vector<int>);
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

int PSICluster::getSize() {
	// Return this cluster's size
	return size;
}

void PSICluster::setReactionNetwork(
		const std::shared_ptr<ReactionNetwork> reactionNetwork) {

	// Call the superclass's method
	Reactant::setReactionNetwork(reactionNetwork);

	// Extract properties from the network
	auto properties = network->getProperties();
	int connectivityLength = network->getAll()->size();

	// Get the enabled reaction type flags
	bool reactionsEnabled = ((*properties)["reactionsEnabled"] == "true");
	bool dissociationsEnabled =
			((*properties)["dissociationsEnabled"] == "true");

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

double PSICluster::getDissociationFlux(const double temperature) {

	int nReactants = 0, oneIndex = -1;
	double diss = 0.0, conc = 0.0;
	std::shared_ptr<Reactant> first, second;
	std::shared_ptr<std::vector<std::shared_ptr<Reactant>> > reactants;
	std::map<std::string, int> oneHe, oneV, oneI;

	// Only try this if the network is available
	if (network != NULL) {
		// Get the total number of reactants in the network
		reactants = network->reactants;
		nReactants = reactants->size();
		// Set the cluster map data for 1 of each species
		oneHe["He"] = 1;
		oneHe["V"] = 0;
		oneHe["I"] = 0;
		oneV["He"] = 0;
		oneV["V"] = 1;
		oneV["I"] = 0;
		oneI["He"] = 0;
		oneI["V"] = 0;
		oneI["I"] = 1;

		// Get this PSICluster or subclasses' cluster map
		std::map<std::string, int> thisMap = getClusterMap();

		// Get the number of species to determine if this
		// cluster is mixed or single
		int numSpecies = (thisMap["He"] > 0) + (thisMap["V"] > 0)
				+ (thisMap["I"] > 0);

		// If no species, throw error
		if (numSpecies == 0) {
			// Bad if we have no species
			throw std::string("Cluster map contains no species");

		} else if (numSpecies == 1) {

			// We know we are a single species,
			// but we need to know which one so we can
			// get the correct species He, V, or I to calculate
			// the dissociation constant.
			if (thisMap["He"]) {
				oneIndex = network->toClusterIndex(oneHe);
			} else if (thisMap["V"]) {
				oneIndex = network->toClusterIndex(oneV);
			} else if (thisMap["I"]) {
				oneIndex = network->toClusterIndex(oneI);
			}
			second = reactants->at(oneIndex);
			// Loop over all reactants and see if we
			// have a dissociation connection
			for (int i = 0; i < nReactants; i++) {
				// Only calculate if we are connected
				if (dissociationConnectivity.at(i) == 1) {
					first = reactants->at(i);
					// Calculate the dissociation flux
					diss += calculateDissociationConstant(first, second,
							temperature)
							* network->reactants->at(i)->getConcentration();
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

double PSICluster::getProductionFlux(const double temperature) {

	// Local declarations
	double flux = 0.0, kPlus = 0.0;
	int thisClusterIndex = 0;
	std::shared_ptr<Reactant> firstReactant, secondReactant;
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
			flux += calculateReactionRateConstant(firstReactant, secondReactant,
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
double PSICluster::getCombinationFlux(const double temperature) {
	// Local declarations
	double flux = 0.0, kPlus = 0.0;
	int thisClusterIndex = 0;
	std::shared_ptr<Reactant> outerReactant, thisReactant;
	int nReactants = 0;

	// Only try this if the network is available
	if (network != NULL) {
		// Set the total network nReactants
		nReactants = combiningReactants.size();
		// Get the index of this cluster in the network
		thisClusterIndex = network->toClusterIndex(getClusterMap());
		thisReactant = network->reactants->at(thisClusterIndex);
		// Loop over all possible clusters
		for (int j = 0; j < nReactants; j++) {
			outerReactant = combiningReactants.at(j);
			// Calculate Second term of production flux
			flux += calculateReactionRateConstant(thisReactant, outerReactant,
					temperature) * outerReactant->getConcentration();
		}

	}

	// Return the production flux
	return (flux * getConcentration());
}

double PSICluster::getTotalFlux(const double temperature) {
	return getProductionFlux(temperature) - getCombinationFlux(temperature)
			+ getDissociationFlux(temperature);
}

double PSICluster::getDiffusionFactor() {
	// Return the diffusion factor
	return diffusionFactor;
}

void PSICluster::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	return;
}

double PSICluster::getDiffusionCoefficient(const double temperature) {
	// Use the Arrhenius equation to compute the diffusion coefficient
	double k_b = xolotlCore::kBoltzmann;
	double kernel = -migrationEnergy / (k_b * temperature);
	return diffusionFactor * exp(kernel);
}

std::vector<double> PSICluster::getBindingEnergies() {
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

double PSICluster::getMigrationEnergy() {
	// Return the migration energy
	return migrationEnergy;
}

void PSICluster::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;
	return;
}

double PSICluster::calculateReactionRateConstant(
		const std::shared_ptr<xolotlCore::Reactant> & firstReactant,
		const std::shared_ptr<xolotlCore::Reactant> & secondReactant,
		const double temperature) {

	// Get the reaction radii
	double r_first =
			(std::dynamic_pointer_cast<PSICluster>(firstReactant))->getReactionRadius();
	double r_second =
			(std::dynamic_pointer_cast<PSICluster>(secondReactant))->getReactionRadius();

	// Get the diffusion coefficients
	double firstDiffusion = (std::dynamic_pointer_cast<PSICluster>(
			firstReactant))->getDiffusionCoefficient(temperature);
	double secondDiffusion = (std::dynamic_pointer_cast<PSICluster>(
			secondReactant))->getDiffusionCoefficient(temperature);

	// Calculate and return
	double k_plus = 4.0 * xolotlCore::pi * (r_first + r_second)
			* (firstDiffusion + secondDiffusion);
	return k_plus;
}

double PSICluster::calculateDissociationConstant(
		const std::shared_ptr<xolotlCore::Reactant> & firstReactant,
		const std::shared_ptr<xolotlCore::Reactant> & secondReactant,
		const double temperature) {

	// Local Declarations
	int bindingEnergyIndex = -1;
	double ra = 1, rb = 1;
	double atomicVolume = 1.0;
	std::shared_ptr<PSICluster> castedSecondReactant =
			(std::dynamic_pointer_cast<PSICluster>(secondReactant));
	std::map<std::string, int> clusterMap =
			castedSecondReactant->getClusterMap();

	// Get the binding energy index
	if (clusterMap["He"] == 1 && clusterMap["V"] == 0 && clusterMap["I"] == 0) {
		bindingEnergyIndex = 0;
	} else if (clusterMap["He"] == 0 && clusterMap["V"] == 1
			&& clusterMap["I"] == 0) {
		bindingEnergyIndex = 1;
	} else if (clusterMap["He"] == 0 && clusterMap["V"] == 0
			&& clusterMap["I"] == 1) {
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

bool PSICluster::isProductReactant(const int reactantI, const int reactantJ) {
	// Base class should just return false
	return false;
}

const double PSICluster::getReactionRadius() const {
	return 0.0;
}

std::shared_ptr<std::vector<int>> PSICluster::getConnectivity() {

	int connectivityLength = network->reactants->size();

	// The reaction and dissociate vectors must be the same length
	// as the number of reactants
	if (reactionConnectivity.size() != connectivityLength) {
		throw std::string("The reaction vector is an incorrect length");
	}

	if (dissociationConnectivity.size() != connectivityLength) {
		throw std::string("The dissociation vector is an incorrect length");
	}

	// Resize the array if required
	if (connectivityLength != connectivity->size()) {
		connectivity->resize(connectivityLength, 0);
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
	reactionConnectivity.resize(network->reactants->size(), 0);
}

void PSICluster::createDissociationConnectivity() {
	// By default, generate an array with a zero for each reactant
	// in the network
	dissociationConnectivity.clear();
	dissociationConnectivity.resize(network->reactants->size(), 0);
}

std::map<std::string, int> PSICluster::getClusterMap() {
	// Create an empty cluster map
	std::map<std::string, int> clusterMap;
	return clusterMap;
}
