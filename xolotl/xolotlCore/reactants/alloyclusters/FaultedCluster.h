#ifndef FAULTEDCLUSTER_H
#define FAULTEDCLUSTER_H

// Includes
#include "AlloyCluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial.
 */
class FaultedCluster: public AlloyCluster {

private:

	/**
	 * The default constructor is private because NEClusters must always be
	 * initialized with a size and performance handler registry
	 */
	FaultedCluster() :
		AlloyCluster() {}

public:

	/**
	 * The constructor. All FaultedClusters must be initialized with a size.
	 *
	 * @param n The size of the cluster
	 * @param registry The performance handler registry
	 */
	FaultedCluster(int n, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Destructor
	 */
	~FaultedCluster() {}

	/**
	 * Returns a reactant created using the copy constructor
	 */
	virtual std::shared_ptr<IReactant> clone() {
		return std::shared_ptr<IReactant> (new FaultedCluster(*this));
	}

}; //end class FaultedCluster

} /* namespace xolotlCore */
#endif
