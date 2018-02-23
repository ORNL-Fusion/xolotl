#ifndef NEXECLUSTER_H
#define NEXECLUSTER_H

// Includes
#include "NECluster.h"
#include <xolotlPerf.h>
#include <Constants.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of xenon.
 */
class NEXeCluster: public NECluster {

private:
	static std::string buildName(IReactant::SizeType nXe) {
		std::stringstream nameStream;
		nameStream << "Xe_" << nXe;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	NEXeCluster() = delete;

	/**
	 * The constructor. All NEXeClusters must be initialized with a size.
	 *
	 * @param nXe the number of xenon atoms in the cluster
	 * @param registry The performance handler registry
	 */
	NEXeCluster(int nXe, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			NECluster(_network, registry, buildName(nXe)) {

		// Set the size
		size = nXe;
		// Update the composition map
		composition[toCompIdx(Species::Xe)] = size;

		// Set the typename appropriately
		type = ReactantType::Xe;

		// Compute the reaction radius
		double FourPi = 4.0 * xolotlCore::pi;
		reactionRadius = 1.05
				* pow((3.0 * 85.0 * (double) size) / FourPi, (1.0 / 3.0))
				/ 10.0;
		if (size == 1)
			reactionRadius = 0.3;

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	NEXeCluster(const NEXeCluster& other) = delete;

	/**
	 * Destructor
	 */
	~NEXeCluster() {
	}

};
//end class NEXeCluster

} /* namespace xolotlCore */
#endif
