#ifndef ALLOYINTCLUSTER_H
#define ALLOYINTCLUSTER_H

// Includes
#include "AlloyCluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial.
 */
class AlloyIntCluster: public AlloyCluster {

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	AlloyIntCluster() = delete;

	/**
	 * The constructor. All AlloyIntClusters must be initialized with a size.
	 *
	 * @param n The size of the cluster
	 * @param registry The performance handler registry
	 */
	AlloyIntCluster(int n, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			AlloyCluster(_network, registry) {
		// Set the size
		size = n;
		// Update the composition map
		composition[toCompIdx(Species::I)] = size;

		// Set the reactant name appropriately
		std::stringstream nameStream;
		nameStream << "I_" << size;
		name = nameStream.str();
		// Set the typename appropriately
		type = ReactantType::I;

		// Define the diffusion pre-factor
		{
			double jumpDistance = xolotlCore::alloyLatticeConstant / sqrt(2.0);
			double phononFrequency = 9.6e12;
			double jumpsPerPhonon = 1.0;
			double prefactorExponent = -1.0;
			diffusionFactor = phononFrequency * jumpsPerPhonon * jumpDistance
					* jumpDistance * pow(double(size), prefactorExponent)
					/ (6.0);
		}

		// Define the formation energy
		formationEnergy = _network.getFormationEnergy(type, size);

		// Define the migration energy
		migrationEnergy = 0.5;

		// Define the reaction radius (using spherical)
		reactionRadius = _network.getReactionRadius(type, size);

		return;
	}

	/**
	 * Destructor
	 */
	~AlloyIntCluster() {
	}

	// Int clusters are considered spheres so set isSphere to true
	bool isSphere() const override {
		return true;
	}

	double getEmissionFlux(int i) const {
		// Initial declarations
		double flux = AlloyCluster::getEmissionFlux(i);

		// Compute the loss to dislocation sinks
		// Including 20% interstitial bias
		flux += 1.2 * xolotlCore::alloysinkStrength * diffusionCoefficient[i]
			* concentration;

		return flux;
	}

	/**
	* This operation computes the partial derivatives due to emission
	* reactions.
	*
	* @param partials The vector into which the partial derivatives should be
	* inserted. This vector should have a length equal to the size of the
	* network.
	* @param i The location on the grid in the depth direction
	*/
	void getEmissionPartialDerivatives(std::vector<double> & partials,
			int i) const {
		// Initial declarations
		AlloyCluster::getEmissionPartialDerivatives(partials, i);

		// Compute the loss to dislocation sinks
		// k^2 * D * C, including 20% interstitial bias
		partials[id + 1] -= 1.2 * xolotlCore::alloysinkStrength * diffusionCoefficient[i];

		return;
	}

};
//end class AlloyIntCluster

} /* namespace xolotlCore */
#endif
