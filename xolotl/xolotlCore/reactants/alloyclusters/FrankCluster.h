#ifndef FRANKCLUSTER_H
#define FRANKCLUSTER_H

// Includes
#include "AlloyCluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial.
 */
class FrankCluster: public AlloyCluster {

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	FrankCluster() = delete;

	/**
	 * The constructor. All FrankClusters must be initialized with a size.
	 *
	 * @param n The size of the cluster
	 * @param registry The performance handler registry
	 */
	FrankCluster(int n, IReactionNetwork& _network,
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
			//double jumpDistance = xolotlCore::alloyLatticeConstant / sqrt(2.0);
			//double phononFrequency = 1.0e13;
			//double jumpsPerPhonon = 1.0;
			//double prefactorExponent = -1.0;
			//diffusionFactor = phononFrequency * jumpsPerPhonon * jumpDistance
			//    * jumpDistance * pow(double(size),prefactorExponent) / (6.0);
			diffusionFactor = 0.0;
		}

		// Define the formation energy
		{
			// USING PLACEHOLDER VALUES => UPDATE WITH CORRECT VALUES
			double A = 4.0;
			double B = 2.0;
			formationEnergy = A + B * (pow(double(size), 2.0 / 3.0) - 1.0);
		}

		// Define the migration energy
		migrationEnergy = 0.5;

		// Define the reaction radius
		reactionRadius = 0.5 * xolotlCore::alloyLatticeConstant
				* sqrt(double(size) * sqrt(3.0) / xolotlCore::pi);
	}

	/**
	 * Destructor
	 */
	~FrankCluster() {
	}

};
//end class FrankCluster

} /* namespace xolotlCore */
#endif
