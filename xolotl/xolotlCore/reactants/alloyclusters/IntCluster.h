#ifndef INTCLUSTER_H
#define INTCLUSTER_H

// Includes
#include "AlloyCluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial.
 */
class IntCluster: public AlloyCluster {

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	IntCluster() = delete;

	/**
	 * The constructor. All IntClusters must be initialized with a size.
	 *
	 * @param n The size of the cluster
	 * @param registry The performance handler registry
	 */
	IntCluster(int n, IReactionNetwork& _network,
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
			double phononFrequency = 1.0e13;
			double jumpsPerPhonon = 1.0;
			double prefactorExponent = -1.0;
			diffusionFactor = phononFrequency * jumpsPerPhonon * jumpDistance
					* jumpDistance * pow(double(size), prefactorExponent)
					/ (6.0);
		}

		// Define the formation energy
		{
			// USING PLACEHOLDER VALUES => UPDATE WITH CORRECT VALUES
			double A = 4.0;
			double B = 3.5;
			formationEnergy = A + B * (pow(double(size), 2.0 / 3.0) - 1.0);
		}

		// Define the migration energy
		migrationEnergy = 0.5;

		// Define the reaction radius (using spherical)
		reactionRadius = 0.5 * xolotlCore::alloyLatticeConstant
				* pow(1.5 * double(size) / xolotlCore::pi, 1.0 / 3.0);

		return;
	}

	/**
	 * Destructor
	 */
	~IntCluster() {
	}

};
//end class IntCluster

} /* namespace xolotlCore */
#endif
