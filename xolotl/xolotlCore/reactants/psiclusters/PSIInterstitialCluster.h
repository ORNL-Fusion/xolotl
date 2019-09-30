#ifndef PSIINTERSTITIALCLUSTER_H
#define PSIINTERSTITIALCLUSTER_H

// Includes
#include <sstream>
#include "PSICluster.h"
#include <Constants.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial defects.
 */
class PSIInterstitialCluster: public PSICluster {

	static std::string buildName(IReactant::SizeType nI) {
		// Set the reactant name appropriately
		std::stringstream nameStream;
		nameStream << "I_" << nI;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSIInterstitialCluster() = delete;

	/**
	 * The constructor. All PSIInterstitialClusters must be initialized with
	 * a size.
	 *
	 * @param nI The number of interstitial defect in this cluster
	 * @param registry The performance handler registry
	 */
	PSIInterstitialCluster(int nI, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSICluster(_network, registry, buildName(nI)) {
		// Set the size
		size = nI;
		// Update the composition map
		composition[toCompIdx(Species::I)] = size;

		// Set the typename appropriately
		type = ReactantType::I;

		// Compute the reaction radius
		double EightPi = 8.0 * xolotlCore::pi;
		double aCubed = pow(xolotlCore::tungstenLatticeConstant, 3.0);
		double termOne = 1.15 * (sqrt(3.0) / 4.0)
				* xolotlCore::tungstenLatticeConstant;
		double termTwo = pow((3.0 / EightPi) * aCubed * size, (1.0 / 3.0));
		double termThree = pow((3.0 / EightPi) * aCubed, (1.0 / 3.0));
		reactionRadius = termOne + termTwo - termThree;

		// Bounds on He, D, T, and V
		bounds[0] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		bounds[1] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		bounds[2] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		bounds[3] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	PSIInterstitialCluster(const PSIInterstitialCluster& other) = delete;

	/**
	 * The Destructor
	 */
	~PSIInterstitialCluster() {
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
//end class PSIInterstitialCluster

} /* end namespace xolotlCore */
#endif
