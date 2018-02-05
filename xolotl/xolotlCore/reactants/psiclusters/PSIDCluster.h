#ifndef PSIDCLUSTER_H
#define PSIDCLUSTER_H

// Includes
#include "PSIHCluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of deuterium.
class PSIDCluster: public PSIHCluster {

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSIDCluster() = delete;

	//! The constructor. All PSIDClusters must be initialized with a size.
	PSIDCluster(int nD, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSIHCluster(nD, _network, registry, "Deuterium") {
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	PSIDCluster(const PSIDCluster& other) = delete;

	//! Destructor
	~PSIDCluster() {
	}

};
//end class PSIDCluster

} /* end namespace xolotlCore */
#endif
