// Includes
#include "SuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

SuperCluster::SuperCluster(double numHe, double numV, int width, double radius,
		double energy,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(1, registry), numHe(numHe), numV(numV) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = (int) (numHe + numV);

	// Update the composition map
	compositionMap[heType] = (int) (numHe * (double) width);
	compositionMap[vType] = (int) (numV * (double) width);

	// Set the width
	sectionWidth = width;

	// Set the reaction radius and formation energy
	reactionRadius = radius;
	formationEnergy = energy;

	// Set the diffusion factor and the migration energy
	migrationEnergy = std::numeric_limits<double>::infinity();
	diffusionFactor = 0.0;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << numHe << "V_" << numV;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "Super";

	return;
}

SuperCluster::SuperCluster(const SuperCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numV = other.numV;

	return;
}

std::shared_ptr<Reactant> SuperCluster::clone() {
	std::shared_ptr<Reactant> reactant(new SuperCluster(*this));

	return reactant;
}

void SuperCluster::createReactionConnectivity() {
	// Aggregate the reacting pairs and combining reactants from the heVVector
	// Loop on the heVVector
	for (int i = 0; i < heVVector.size(); i++) {
		// Get both production vectors
		auto react = heVVector[i]->reactingPairs;
		auto combi = heVVector[i]->combiningReactants;

		// Loop on its reacting pairs
		for (int l = 0; l < react.size(); l++) {
			// Add the pairs to our vector
			reactingPairs.push_back(react[l]);
		}

		// Loop on its combining reactants
		for (int l = 0; l < combi.size(); l++) {
			// Add the pairs to our vector
			combiningReactants.push_back(combi[l]);
		}
	}

	return;
}

void SuperCluster::createDissociationConnectivity() {
	// Aggregate the dissociating and emission pairs from the heVVector
	// Loop on the heVVector
	for (int i = 0; i < heVVector.size(); i++) {
		// Get both dissociation vectors
		auto disso = heVVector[i]->dissociatingPairs;
		auto emi = heVVector[i]->emissionPairs;

		// Loop on its dissociating pairs
		for (int l = 0; l < disso.size(); l++) {
			// Add the pairs to our vector
			dissociatingPairs.push_back(disso[l]);
		}

		// Loop on its emission pairs
		for (int l = 0; l < emi.size(); l++) {
			// Add the pairs to our vector
			emissionPairs.push_back(emi[l]);
		}
	}

	return;
}

void SuperCluster::computeRateConstants() {
	// Initialize all the effective vectors
	effReactingPairs.clear();
	effCombiningReactants.clear();
	effDissociatingPairs.clear();
	effEmissionPairs.clear();

	// Local declarations
	PSICluster *firstReactant, *secondReactant, *combiningReactant,
		*dissociatingCluster, *otherEmittedCluster, *firstCluster,
		*secondCluster;
	double rate = 0.0;

	// Compute the reaction constant associated to the reacting pairs
	// Set the total number of reacting pairs
	int nPairs = reactingPairs.size();
	// Loop on them
	for (int i = 0; i < nPairs; i++) {
		// Get the reactants
		firstReactant = reactingPairs[i].first;
		secondReactant = reactingPairs[i].second;
		// Compute the reaction constant
		rate = calculateReactionRateConstant(*firstReactant,
				*secondReactant);
		// Set it in the pair
		reactingPairs[i].kConstant = rate / (double) sectionWidth;

		// Add the reacting pair to the effective vector
		// if the rate is not 0.0
		if (!xolotlCore::equal(rate, 0.0)) {
			effReactingPairs.push_back(&reactingPairs[i]);
		}
	}

	// Compute the reaction constant associated to the combining reactants
	// Set the total number of combining reactants
	int nReactants = combiningReactants.size();
	// Loop on them
	for (int i = 0; i < nReactants; i++) {
		// Get the reactants
		combiningReactant = combiningReactants[i].combining;
		// Compute the reaction constant
		rate = calculateReactionRateConstant(*this, *combiningReactant);
		// Set it in the combining cluster
		combiningReactants[i].kConstant = rate / (double) sectionWidth;

		// Add the combining reactant to the effective vector
		// if the rate is not 0.0
		if (!xolotlCore::equal(rate, 0.0)) {
			effCombiningReactants.push_back(&combiningReactants[i]);

			// Add itself to the list again to account for the correct rate
			if (id == combiningReactant->getId())
				effCombiningReactants.push_back(&combiningReactants[i]);
		}
	}

	// Compute the dissociation constant associated to the dissociating clusters
	// Set the total number of dissociating clusters
	nPairs = dissociatingPairs.size();
	// Loop on them
	for (int i = 0; i < nPairs; i++) {
		dissociatingCluster = dissociatingPairs[i].first;
		// The second element of the pair is the cluster that is also
		// emitted by the dissociation
		otherEmittedCluster = dissociatingPairs[i].second;
		// Compute the dissociation constant
		// The order of the cluster is important here because of the binding
		// energy used in the computation. It is taken from the type of the first cluster
		// which must be the single one
		if (size == 1) {
			// "this" is the single size one
			rate = calculateDissociationConstant(*dissociatingCluster, *this,
					*otherEmittedCluster);
		} else {
			// otherEmittedCluster is the single size one
			rate = calculateDissociationConstant(*dissociatingCluster,
					*otherEmittedCluster, *this);

		}
		// Set it in the pair
		dissociatingPairs[i].kConstant = rate / (double) sectionWidth;

		// Add the dissociating pair to the effective vector
		// if the rate is not 0.0
		if (!xolotlCore::equal(rate, 0.0)) {
			effDissociatingPairs.push_back(&dissociatingPairs[i]);

			// Add itself to the list again to account for the correct rate
			if (id == otherEmittedCluster->getId())
				effDissociatingPairs.push_back(&dissociatingPairs[i]);
		}
	}

	// Compute the dissociation constant associated to the emission of pairs of clusters
	// Set the total number of emission pairs
	nPairs = emissionPairs.size();
	// Loop on them
	for (int i = 0; i < nPairs; i++) {
		firstCluster = emissionPairs[i].first;
		secondCluster = emissionPairs[i].second;
		// Compute the dissociation rate
		rate = calculateDissociationConstant(*this, *firstCluster,
				*secondCluster);
		// Set it in the pair
		emissionPairs[i].kConstant = rate / (double) sectionWidth;

		// Add the emission pair to the effective vector
		// if the rate is not 0.0
		if (!xolotlCore::equal(rate, 0.0)) {
			effEmissionPairs.push_back(&emissionPairs[i]);
		}
	}

	// Shrink the arrays to save some space
	effReactingPairs.shrink_to_fit();
	effCombiningReactants.shrink_to_fit();
	effDissociatingPairs.shrink_to_fit();
	effEmissionPairs.shrink_to_fit();

	return;
}
