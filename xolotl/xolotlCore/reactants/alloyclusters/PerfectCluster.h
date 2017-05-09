#ifndef PERFECTCLUSTER_H
#define PERFECTCLUSTER_H

// Includes
#include "AlloyCluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial.
 */
class PerfectCluster: public AlloyCluster {

private:

	/**
	 * The default constructor is private because NEClusters must always be
	 * initialized with a size and performance handler registry
	 */
	PerfectCluster() :
		AlloyCluster() {}

public:

	/**
	 * The constructor. All PerfectClusters must be initialized with a size.
	 *
	 * @param n The size of the cluster
	 * @param registry The performance handler registry
	 */
	PerfectCluster(int n, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Destructor
	 */
	~PerfectCluster() {}

	/**
	 * Returns a reactant created using the copy constructor
	 */
	virtual std::shared_ptr<IReactant> clone() {
		return std::shared_ptr<IReactant> (new PerfectCluster(*this));
	}


}; //end class PerfectCluster

} /* namespace xolotlCore */
#endif
