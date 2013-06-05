#include "PSICluster.h"
#include <Constants.h>
#include <math.h>

using namespace xolotlCore;

PSICluster::PSICluster(const int clusterSize) : Reactant() {

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

PSICluster::~PSICluster() {}

double PSICluster::getDissociationFlux() {
	return 0.0;
}

double PSICluster::getProductionFlux() {
	return 0.0;
}

int PSICluster::getSize() {
	return size;
}

double PSICluster::getGenByCapt() {
	return 0.0;
}

double PSICluster::getGenByAnn() {
	return 0.0;
}

double PSICluster::getDiffusionFactor() {
	//TODO Auto-generated method stub
	return diffusionFactor;
}

void PSICluster::setDiffusionFactor(const double factor) {

	diffusionFactor = factor;

	return;
}

double PSICluster::getDiffusionCoefficient(const double temperature) {
	// Use the Arrhenius equation to compute the diffusion coefficient
	double k_b = xolotlCore::kBoltzmann;
	double kernel = -migrationEnergy/(k_b*temperature);
	return diffusionFactor*exp(kernel);
}

std::vector<double> PSICluster::getBindingEnergies() {
	// Local Declarations
	std::vector<double> energyVector;

	energyVector = bindingEnergies;

	return energyVector;
}

void PSICluster::setBindingEnergies(const std::vector<double> energies) {
	bindingEnergies = energies;
	return;
}

double PSICluster::getMigrationEnergy() {
	//TODO Auto-generated method stub
	return migrationEnergy;
}

void PSICluster::setMigrationEnergy(const double energy) {

	migrationEnergy = energy;

	return;
}
