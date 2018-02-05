#ifndef PSIHCLUSTER_H
#define PSIHCLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of hydrogen.
class PSIHCluster: public PSICluster {

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSIHCluster() = delete;

	/**
	 * The constructor. All PSIHClusters must be initialized with a size.
	 * @param nH the number of hydrogen atoms in the cluster
	 */
	PSIHCluster(int nH, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry,
			const std::string& _name = "Hydrogen") :
			PSICluster(_network, registry, _name) {

		// Set the size appropriately
		size = nH;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	PSIHCluster(const PSIHCluster& other) = delete;

	//! The destructor
	~PSIHCluster() {
	}

};
//end class PSIHCluster

} /* end namespace xolotlCore */
#endif
