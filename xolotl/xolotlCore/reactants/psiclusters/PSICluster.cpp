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

PSICluster::PSICluster(const PSICluster &other) :
		Reactant(other), size(other.size), diffusionFactor(
				other.diffusionFactor), bindingEnergies(other.bindingEnergies), migrationEnergy(
				other.migrationEnergy) {
}

PSICluster::~PSICluster() {
}

double PSICluster::getTotalFlux(const double temperature) {
	return getProductionFlux(temperature) + getDissociationFlux(temperature);
}

double PSICluster::getDissociationFlux(const double temperature) {

	int nReactants = network->reactants->size(), oneIndex = -1;
	double diss = 0.0, conc = 0.0;
	std::vector<int> dissConnections = getDissociationConnectivity();
	std::map<std::string, int> oneHe, oneV, oneI;

	// Set the cluster map data for 1 of each species
	oneHe["He"] = 1; oneHe["V"] = 0; oneHe["I"] = 0;
	oneV["He"] = 0; oneV["V"] = 1; oneV["I"] = 0;
	oneI["He"] = 0; oneI["V"] = 0; oneI["I"] = 1;

	// Get this PSICluster or subclasses' cluster map
	std::map<std::string, int> thisMap = getClusterMap();

	// Get the number of species to determine if this
	// cluster is mixed or single
	int numSpecies = (thisMap["He"] > 0) + (thisMap["V"] > 0) + (thisMap["I"] > 0);

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

		// Loop over all reactants and see if we
		// have a dissociation connection
		for (int i = 0; i < nReactants; i++) {
			// Only calculate if we are connected
			if (dissConnections.at(i) == 1) {
				// Calculate the dissociation flux
				diss = diss + calculateDissociationConstant(i, oneIndex,
								temperature)
								* network->reactants->at(i)->getConcentration();
			}
		}
	} else if (numSpecies == 2) {
		throw std::string("Mixed Species dissociation flux must be implemented by subclass.");
	}

	// Return the flux
	return diss;
}

double PSICluster::getProductionFlux(const double temperature) {
	// Local declarations
	double fluxOne = 0.0, fluxTwo = 0.0, kPlus = 0.0;
	int thisClusterIndex = 0;

	// This cluster's index in the reactants array - this is Andrew's
	thisClusterIndex = network->toClusterIndex(getClusterMap());

	// Loop over all possible clusters
	for (int j = 0; j < network->reactants->size(); j++) {
		for (int k = 0; k < network->reactants->size(); k++) {
			// If the jth and kth reactants react to produce this reactant...
			if (network->isConnected(j, k) && isProductReactant(j, k)) {
				// This fluxOne term considers all reactions that
				// produce C_i
				fluxOne = fluxOne
						+ calculateReactionRateConstant(j, k, temperature)
								* network->reactants->at(j)->getConcentration()
								* network->reactants->at(k)->getConcentration();
			}
		}

		// Calculate Second term of production flux
		// this acts to take away from the current reactant
		// as it is reacting with others, thus decreasing itself.
		// This considers all populations that are produced by C_i
		if (network->isConnected(j, thisClusterIndex)) {
			fluxTwo = fluxTwo
					+ calculateReactionRateConstant(thisClusterIndex, j,
							temperature)
							* network->reactants->at(j)->getConcentration();
		}
	}

	// Return the production flux
	return fluxOne - (fluxTwo * getConcentration());
}

int PSICluster::getSize() {
	// Return this cluster's size
	return size;
}

double PSICluster::getGenByCapt() {
	return 0.0;
}

double PSICluster::getGenByAnn() {
	return 0.0;
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

double PSICluster::calculateReactionRateConstant(int i, int j, const double temperature) {

	// Get the reaction radii and diffusion coefficients
	double ra = std::dynamic_pointer_cast<PSICluster>(network->reactants->at(i))->getReactionRadius();
	double rb = std::dynamic_pointer_cast<PSICluster>(network->reactants->at(j))->getReactionRadius();

	// Get the Diffusion coefficients
	double iDiffusion = std::dynamic_pointer_cast<PSICluster>(
			network->reactants->at(i))->getDiffusionCoefficient(temperature);
	double jDiffusion = std::dynamic_pointer_cast<PSICluster>(
			network->reactants->at(j))->getDiffusionCoefficient(temperature);

	// Calculate and return
	return 4 * xolotlCore::pi * (ra + rb) * (iDiffusion + jDiffusion);
}

double PSICluster::calculateDissociationConstant(int i, int j, double temperature) {

	// Local Declarations
	int bindingEnergyIndex = -1;
	double ra = 1, rb = 1;
	double atomicVolume = 1.0;
	std::map<std::string, int> clusterMap = network->toClusterMap(j);

	// Get the binding energy index
	if (clusterMap["He"] == 1 && clusterMap["V"] == 0 && clusterMap["I"] == 0) {
		bindingEnergyIndex = 0;
	} else if (clusterMap["He"] == 0 && clusterMap["V"] == 1 && clusterMap["I"] == 0) {
		bindingEnergyIndex = 1;
	} else if (clusterMap["He"] == 0 && clusterMap["V"] == 0 && clusterMap["I"] == 1) {
		bindingEnergyIndex = 2;
	} else {
		return 0.0;
	}

	// Calculate the Reaction Rate Constant -- Cant use this, weird indices change in paper
	double kPlus = calculateReactionRateConstant(i, j, temperature);

	// Calculate and return
	return (1 / atomicVolume) * kPlus
			* exp(bindingEnergies.at(bindingEnergyIndex)
							/ (xolotlCore::kBoltzmann * temperature));
}

bool PSICluster::isProductReactant(int reactantI, int reactantJ) {
	// Base class should just return false
	return false;
}

double PSICluster::getReactionRadius() {
	return 0.0;
}


std::vector<int> PSICluster::getConnectivity() {
	
	std::vector<int> reactionConn = getReactionConnectivity();
	std::vector<int> dissConn = getDissociationConnectivity();
	
	// The vectors must be the same length.
	
	if (reactionConn.size() != dissConn.size()) {
		throw std::string("The reaction and dissociation vectors "
			"must be the same length");
	}
	
	// Merge the two vectors such that the final vector contains
	// a 1 at a positioin if either of the connectivity arrays
	// have a 1
	
	int connLength = reactionConn.size();
	for (int i = 0; i < connLength; i++) {
		
		// We can modify the reaction vector in place since it
		// will not be needed by this method again.
		reactionConn[i] |= dissConn[i];
	}
	
	return reactionConn;
}


std::vector<int> PSICluster::getReactionConnectivity() {
	// By default, return an array with a zero for each reactant
	// in the network
	
	std::vector<int> reactionConn(network->reactants->size(), 0);
	return reactionConn;
}


std::vector<int> PSICluster::getDissociationConnectivity() {
	// By default, return an array with a zero for each reactant
	// in the network
	
	std::vector<int> dissConn(network->reactants->size(), 0);
	return dissConn;
}

std::map<std::string, int> PSICluster::getClusterMap() {
	std::map<std::string, int> dummy;
	return dummy;
}
