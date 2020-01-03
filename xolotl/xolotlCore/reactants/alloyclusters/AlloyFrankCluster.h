#ifndef ALLOYFRANKCLUSTER_H
#define ALLOYFRANKCLUSTER_H

// Includes
#include "AlloyCluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial.
 */
class AlloyFrankCluster: public AlloyCluster {

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	AlloyFrankCluster() = delete;

	/**
	 * The constructor. All AlloyFrankClusters must be initialized with a size.
	 *
	 * @param n The size of the cluster
	 * @param registry The performance handler registry
	 */
	AlloyFrankCluster(int n, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			AlloyCluster(_network, registry) {
		// Set the size
		size = n;
		// Update the composition map
		composition[toCompIdx(Species::Frank)] = size;

		// Set the reactant name appropriately
		std::stringstream nameStream;
		nameStream << "Fr_" << size;
		name = nameStream.str();
		// Set the typename appropriately
		type = ReactantType::Frank;

		// Define the diffusion pre-factor
		{
			//double jumpDistance = network.getLatticeParameter() / sqrt(2.0);
			//double phononFrequency = 1.0e13;
			//double jumpsPerPhonon = 1.0;
			//double prefactorExponent = -1.0;
			//diffusionFactor = phononFrequency * jumpsPerPhonon * jumpDistance
			//    * jumpDistance * pow(double(size),prefactorExponent) / (6.0);
			diffusionFactor = 0.0;
		}

		// Define the formation energy
		formationEnergy = _network.getFormationEnergy(type, size);

		// Define the migration energy
		migrationEnergy = 0.5;

		// Define the reaction radius
		reactionRadius = _network.getReactionRadius(type, size);
	}

	/**
	 * Destructor
	 */
	~AlloyFrankCluster() {
	}

	/**
	 * Add grid points to the vector of diffusion coefficients or remove
	 * them if the value is negative.
	 *
	 * @param i The number of grid point to add or remove
	 */
	void addGridPoints(int i) override {
		if (diffusionFactor > 0.0) {
			Reactant::addGridPoints(i);
		}

		// Don't do anything
		return;
	}

	/**
	 * This operation sets the temperature at which the reactant currently
	 * exists. Temperature-dependent quantities are recomputed when this
	 * operation is called, so the temperature should always be set first.
	 *
	 * @param temp The new cluster temperature
	 * @param i The location on the grid
	 */
	void setTemperature(double temp, int i) override{
		if (diffusionFactor > 0.0) {
			Reactant::setTemperature(temp, i);
		}

		// Don't do anything
		return;
	}

	/**
	 * This operation returns the diffusion coefficient for this reactant and is
	 * calculated from the diffusion factor.
	 *
	 * @param i The position on the grid
	 * @return The diffusion coefficient
	 */
	double getDiffusionCoefficient(int i) const override {
		if (diffusionFactor > 0.0) {
			return Reactant::getDiffusionCoefficient(i);
		}

		return 0.0;
	}

};
//end class AlloyFrankCluster

} /* namespace xolotlCore */
#endif
