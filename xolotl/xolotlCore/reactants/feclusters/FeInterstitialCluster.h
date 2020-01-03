#ifndef FEINTERSTITIALCLUSTER_H
#define FEINTERSTITIALCLUSTER_H

// Includes
#include <sstream>
#include "FeCluster.h"
#include <Constants.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial defects.
 */
class FeInterstitialCluster: public FeCluster {

	static std::string buildName(IReactant::SizeType nI) {
		// Set the reactant name appropriately
		std::stringstream nameStream;
		nameStream << "I_" << nI;
		return nameStream.str();
	}

	/**
	 * The sink strength
	 */
	double sinkStrength;

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	FeInterstitialCluster() = delete;

	/**
	 * The constructor. All FeInterstitialClusters must be initialized with
	 * a size.
	 *
	 * @param nI The number of interstitial defect in this cluster
	 * @param registry The performance handler registry
	 */
	FeInterstitialCluster(int nI, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			FeCluster(_network, registry, buildName(nI)) {

		// Set the size
		size = nI;
		// Update the composition map
		composition[toCompIdx(Species::I)] = size;

		// Set the typename appropriately
		type = ReactantType::I;

		// Compute the reaction radius
		double latticeParam = network.getLatticeParameter();
		double EightPi = 8.0 * xolotlCore::pi;
		reactionRadius = latticeParam
				* pow((3.0 / EightPi) * size, (1.0 / 3.0));

		// Bounds on He and V
		heBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		vBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));

		//! Parameters for biased sink in the iron case
		double r0 = latticeParam * 0.75 * sqrt(3.0);
		double reactionRadius = latticeParam * cbrt(3.0 / xolotlCore::pi) * 0.5;
		constexpr double rho = 0.0003;
		sinkStrength = 1.05 * -4.0 * xolotlCore::pi * rho
				/ log(
						xolotlCore::pi * rho * (reactionRadius + r0)
								* (reactionRadius + r0));

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	FeInterstitialCluster(const FeInterstitialCluster& other) = delete;

	/**
	 * The Destructor
	 */
	~FeInterstitialCluster() {
	}

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux(int i) override {
		// Initial declarations
		double flux = FeCluster::getEmissionFlux(i);

		// Compute the loss to dislocation sinks
		if (size < 2) {
			// bias * k^2 * D * C
			flux += sinkStrength * diffusionCoefficient[i] * concentration;
		}

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
			int i) const override {
		// Initial declarations
		FeCluster::getEmissionPartialDerivatives(partials, i);

		// Compute the loss to dislocation sinks
		if (size < 2) {
			// bias * k^2 * D * C
			partials[id - 1] -= sinkStrength * diffusionCoefficient[i];
		}

		return;
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
	void setTemperature(double temp, int i) override {
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
//end class FeInterstitialCluster

} /* end namespace xolotlCore */
#endif
