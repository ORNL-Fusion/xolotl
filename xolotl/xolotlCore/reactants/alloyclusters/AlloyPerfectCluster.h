#ifndef ALLOYPERFECTCLUSTER_H
#define ALLOYPERFECTCLUSTER_H

// Includes
#include "AlloyCluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial.
 */
class AlloyPerfectCluster: public AlloyCluster {

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	AlloyPerfectCluster() = delete;

	/**
	 * The constructor. All AlloyPerfectClusters must be initialized with a size.
	 *
	 * @param n The size of the cluster
	 * @param registry The performance handler registry
	 */
	AlloyPerfectCluster(int n, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			AlloyCluster(_network, registry) {
		// Set the size
		size = n;
		// Update the composition map
		composition[toCompIdx(Species::Perfect)] = size;

		// Set the reactant name appropriately
		std::stringstream nameStream;
		nameStream << "Pe_" << size;
		name = nameStream.str();
		// Set the typename appropriately
		type = ReactantType::Perfect;

		// Define the diffusion pre-factor
		if (size < 70) {
			double jumpDistance = xolotlCore::alloyLatticeConstant / sqrt(2.0);
			double phononFrequency = 9.6e12;
			double jumpsPerPhonon = 1.0;
			double prefactorExponent = -1.0;
			diffusionFactor = phononFrequency * jumpsPerPhonon * jumpDistance
					* jumpDistance * pow(double(size), prefactorExponent)
					/ (6.0);
		} else {
			diffusionFactor = 0.0;
		}

		// Define the formation energy
		formationEnergy = _network.getFormationEnergy(type, size);

		// Define the migration energy
		migrationEnergy = 0.7;

		// Define the reaction radius
		reactionRadius = _network.getReactionRadius(type, size);

		return;
	}

	/**
	 * Destructor
	 */
	~AlloyPerfectCluster() {
	}

        double getEmissionFlux(int i) const {
                // Initial declarations
                double flux = AlloyCluster::getEmissionFlux(i);

                // Compute the loss to dislocation sinks
                flux += xolotlCore::alloysinkStrength * diffusionCoefficient[i]      
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
                // k^2 * D * C
                partials[id + 1] -= xolotlCore::alloysinkStrength * diffusionCoefficient[i];

                return;
        }

};
//end class AlloyPerfectCluster

} /* namespace xolotlCore */
#endif
