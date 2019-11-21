#ifndef PSIDCLUSTER_H
#define PSIDCLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of deuterium.
class PSIDCluster: public PSICluster {

private:
	static std::string buildName(IReactant::SizeType size) {
		// Set the reactant name appropriately
		std::stringstream nameStream;
		nameStream << "D_" << size;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSIDCluster() = delete;

	/**
	 * The constructor. All PSIDClusters must be initialized with a size.
	 *
	 * @param nD the number of deuterium atoms in the cluster
	 * @param factor The ratio between the He and H radius
	 * @param registry The performance handler registry
	 */
	PSIDCluster(int nD, double factor, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSICluster(_network, registry, buildName(nD)) {
		// Set the size
		size = nD;
		// Update the composition map
		composition[toCompIdx(Species::D)] = size;
		// Set the typename appropriately
		type = ReactantType::D;

		// Compute the reaction radius
		double FourPi = 4.0 * xolotlCore::pi;
		double aCubed = pow(network.getLatticeParameter(), 3);
		double termOne = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * size,
				(1.0 / 3.0));
		double termTwo = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed,
				(1.0 / 3.0));
		reactionRadius = (0.3 + termOne - termTwo) * factor;

		// Bounds on He, D, T, and V
		bounds[0] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		bounds[1] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(size),
				static_cast<IReactant::SizeType>(size + 1));
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
	PSIDCluster(const PSIDCluster& other) = delete;

	/**
	 * Destructor
	 */
	~PSIDCluster() {
	}

};
//end class PSIDCluster

} /* end namespace xolotlCore */
#endif
