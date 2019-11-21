#ifndef PSITCLUSTER_H
#define PSITCLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of tritium.
class PSITCluster: public PSICluster {

private:
	static std::string buildName(IReactant::SizeType size) {
		// Set the reactant name appropriately
		std::stringstream nameStream;
		nameStream << "T_" << size;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSITCluster() = delete;

	/**
	 * The constructor. All PSITClusters must be initialized with a size.
	 *
	 * @param nT the number of tritium atoms in the cluster
	 * @param factor The ratio between the He and H radius
	 * @param registry The performance handler registry
	 */
	PSITCluster(int nT, double factor, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSICluster(_network, registry, buildName(nT)) {
		// Set the size
		size = nT;
		// Update the composition map
		composition[toCompIdx(Species::T)] = size;
		// Set the typename appropriately
		type = ReactantType::T;

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
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		bounds[2] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(size),
				static_cast<IReactant::SizeType>(size + 1));
		bounds[3] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	PSITCluster(const PSITCluster& other) = delete;

	/**
	 * Destructor
	 */
	~PSITCluster() {
	}

};
//end class PSITCluster

} /* end namespace xolotlCore */
#endif
