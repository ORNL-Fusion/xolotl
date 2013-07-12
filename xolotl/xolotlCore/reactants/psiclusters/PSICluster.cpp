#include "PSICluster.h"
#include <Constants.h>
#include <math.h>
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
	return 0.0;
}

double PSICluster::getProductionFlux(const double temperature) {
	// Local declarations
	double fluxOne = 0.0, fluxTwo = 0.0, kPlus = 0.0;
	int thisClusterIndex = 0;
	int numHeClusters = std::stoi(network->properties->at("numHeClusters"));
	int numVClusters = std::stoi(network->properties->at("numVClusters"));
	int numIClusters = std::stoi(network->properties->at("numIClusters"));
	int numSingleSpeciesClusters = numHeClusters + numVClusters + numIClusters;

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
	double ra = 1, rb = 1; // Will set these later
	return 4 * xolotlCore::pi * (ra + rb)
			* (std::dynamic_pointer_cast<PSICluster>(network->reactants->at(i))->getDiffusionCoefficient(
					temperature)
					+ std::dynamic_pointer_cast<PSICluster>(
							network->reactants->at(j))->getDiffusionCoefficient(
							temperature));
}

double PSICluster::calculateDissociationConstant(int i, int j, double temperature) {

	// Local Declarations
	int bindingEnergyIndex = -1;
	double atomicVolume = 1.0;
	std::map<std::string, int> clusterMap = network->toClusterMap(i);

	// Calculate the Reaction Rate Constant
	double kPlus = calculateReactionRateConstant(i, j, temperature);

	// Get the species at index i so we
	// can get the binding energy index


	// Make sure we found a valid binding energy
	if (bindingEnergyIndex == -1) {
		return 0.0;
	} else {
		return (1 / atomicVolume) * kPlus* exp(
					bindingEnergies.at(bindingEnergyIndex)
							/ (xolotlCore::kBoltzmann * temperature));
	}
}

bool PSICluster::isProductReactant(int reactantI, int reactantJ) {
	// Base class should just return false
	return false;
}

std::map<std::string, int> PSICluster::getClusterMap() {
	std::map<std::string, int> map;
	return map;
}
