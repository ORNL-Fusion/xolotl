#ifndef PSITCLUSTER_H
#define PSITCLUSTER_H

// Includes
#include "PSIHCluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of tritium.
class PSITCluster: public PSIHCluster {

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSITCluster() = delete;

	//! The constructor. All PSITClusters must be initialized with a size.
	PSITCluster(int nT, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSIHCluster(nT, _network, registry, "Tritium") {
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	PSITCluster(const PSITCluster& other) = delete;

	//! Destructor
	~PSITCluster() {
	}

};
//end class PSITCluster

} /* end namespace xolotlCore */
#endif
