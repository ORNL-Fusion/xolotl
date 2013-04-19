#include "PSICluster.h"

using namespace xolotlCore;

PSICluster::PSICluster(int clusterSize) : Reactant() {

	// Set the size
	size = (clusterSize > 0) ? clusterSize : 1;
	// Zero out the binding energies
	bindingEnergies[0] = 0.0;
	bindingEnergies[1] = 0.0;
	bindingEnergies[2] = 0.0;
	bindingEnergies[3] = 0.0;
	// Zero out the diffusion factor
	diffusionFactor = 0.0;
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

double PSICluster::getDiffusionCoefficient() {
	return 0.0;
}

double PSICluster::getDiffusionFactor() {
	//TODO Auto-generated method stub
	return diffusionFactor;
}

void PSICluster::setDiffusionFactor(double factor) {

	diffusionFactor = factor;

	return;
}

std::vector<double> PSICluster::getBindingEnergy() {
	// Local Declarations
	std::vector<double> energyVector;

	energyVector.push_back(bindingEnergies[0]);
	energyVector.push_back(bindingEnergies[1]);
	energyVector.push_back(bindingEnergies[2]);
	energyVector.push_back(bindingEnergies[3]);

	return energyVector;
}

void PSICluster::setBindingEnergy(std::vector<double> energies) {
	//TODO Auto-generated method stub
	return;
}
