#ifndef PSIHECLUSTER_H
#define PSIHECLUSTER_H

// Includes
#include <sstream>
#include "PSICluster.h"
#include <xolotlPerf/xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of helium.
 */
class PSIHeCluster: public PSICluster {

private:
	static std::string buildName(IReactant::SizeType size) {
		// Set the reactant name appropriately
		std::stringstream nameStream;
		nameStream << "He_" << size;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSIHeCluster() = delete;

	/**
	 * The constructor. All PSIHeClusters must be initialized with a size.
	 *
	 * @param nHe the number of helium atoms in the cluster
	 * @param registry The performance handler registry
	 */
	PSIHeCluster(int nHe, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSICluster(_network, registry, buildName(nHe)) {
		// Set the size
		size = nHe;
		// Update the composition map
		composition[toCompIdx(Species::He)] = size;
		// Set the typename appropriately
		type = ReactantType::He;

		// Compute the reaction radius
		double FourPi = 4.0 * xolotlCore::pi;
		double aCubed = pow(network.getLatticeParameter(), 3);
		double termOne = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * size,
				(1.0 / 3.0));
		double termTwo = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed,
				(1.0 / 3.0));
		reactionRadius = network.getImpurityRadius() + termOne - termTwo;

		// Bounds on He, D, T, and V
		bounds[0] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(size),
				static_cast<IReactant::SizeType>(size + 1));
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
	PSIHeCluster(const PSIHeCluster& other) = delete;

	/**
	 * Destructor
	 */
	~PSIHeCluster() {
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
//end class PSIHeCluster

} /* namespace xolotlCore */
#endif
