#ifndef PSIHCLUSTER_H
#define PSIHCLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of hydrogen.
class PSIHCluster: public PSICluster {

private:
	static std::string buildName(IReactant::SizeType size) {
		// Set the reactant name appropriately
		std::stringstream nameStream;
		nameStream << "H_" << size;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSIHCluster() = delete;

	/**
	 * The constructor. All PSIHClusters must be initialized with a size.
	 *
	 * @param nH the number of hydrogen atoms in the cluster
	 * @param registry The performance handler registry
	 */
	PSIHCluster(int nH, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSICluster(_network, registry, buildName(nH)) {
		// Don't do anything because the daughter classes will be in charge
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
	PSIHCluster(const PSIHCluster& other) = delete;

	/**
	 * Destructor
	 */
	~PSIHCluster() {
	}

};
//end class PSIHCluster

} /* end namespace xolotlCore */
#endif
