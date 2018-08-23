#ifndef VACCLUSTER_H
#define VACCLUSTER_H

// Includes
#include "AlloyCluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial.
 */
class VacCluster: public AlloyCluster {

private:

	/**
	 * The default constructor is private because NEClusters must always be
	 * initialized with a size and performance handler registry
	 */
	VacCluster() :
		AlloyCluster() {}

public:

	/**
	 * The constructor. All VacClusters must be initialized with a size.
	 *
	 * @param n The size of the cluster
	 * @param registry The performance handler registry
	 */
	VacCluster(int n, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Destructor
	 */
	~VacCluster() {}

	/**
	 * Returns a reactant created using the copy constructor
	 */
	virtual std::shared_ptr<IReactant> clone() {
		return std::shared_ptr<IReactant> (new VacCluster(*this));
	}

}; //end class VacCluster

} /* namespace xolotlCore */
#endif
