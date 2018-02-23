#ifndef PSITCLUSTER_H
#define PSITCLUSTER_H

// Includes
#include "PSIHCluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of tritium.
class PSITCluster: public PSIHCluster {

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
	 * @param registry The performance handler registry
	 */
	PSITCluster(int nT, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSIHCluster(nT, _network, registry) {
		// Set the size
		size = nT;
		// Update the composition map
		composition[toCompIdx(Species::T)] = size;
		// Set the typename appropriately
		type = ReactantType::T;

		// Compute the reaction radius
		reactionRadius = 0.3;

		// Bounds on He and V
		heBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		vBounds = IntegerRange<IReactant::SizeType>(
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
