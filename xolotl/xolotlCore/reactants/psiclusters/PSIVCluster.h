#ifndef PSIVCLUSTER_H
#define PSIVCLUSTER_H

// Includes
#include "PSICluster.h"
#include <xolotlPerf.h>
#include <Constants.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of atomic vacancies.
 */
class PSIVCluster: public PSICluster {

private:
	static std::string buildName(IReactant::SizeType nV) {
		std::stringstream nameStream;
		nameStream << "V_" << nV;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSIVCluster() = delete;

	/**
	 * The constructor. All PSIVClusters must be initialized with a size.
	 *
	 * @param nV the number of atomic vacancies in the cluster
	 * @param registry The performance handler registry
	 */
	PSIVCluster(int nV, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSICluster(_network, registry, buildName(nV)) {
		// Set the size
		size = nV;
		// Update the composition map
		composition[toCompIdx(Species::V)] = size;
		// Set the typename appropriately
		type = ReactantType::V;

		// Compute the reaction radius
		// It is the same formula for HeV clusters
		reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::tungstenLatticeConstant
				+ pow(
						(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0)
								* size) / (8.0 * xolotlCore::pi), (1.0 / 3.0))
				- pow(
						(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0))
								/ (8.0 * xolotlCore::pi), (1.0 / 3.0));

		// Bounds on He and V
		heBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		vBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(size),
				static_cast<IReactant::SizeType>(size + 1));

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	PSIVCluster(const PSIVCluster& other) = delete;

	//! Destructor
	~PSIVCluster() {
	}

};
//end class PSIVCluster

} /* end namespace xolotlCore */

#endif
