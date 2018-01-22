#ifndef VCLUSTER_H
#define VCLUSTER_H

// Includes
#include "PSICluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of atomic vacancies.
 */
class VCluster: public PSICluster {

private:

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	VCluster() :
		PSICluster()
	{ }

public:

	/**
	 * The constructor. All VClusters must be initialized with a size.
	 *
	 * @param nV the number of atomic vacancies in the cluster
	 * @param registry The performance handler registry
	 */
	VCluster(int nV, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	//! Destructor
	~VCluster() {}

	/**
	 * Returns a reactant created using the copy constructor
	 */
	virtual std::shared_ptr<IReactant> clone() {
		return std::shared_ptr<IReactant> (new VCluster(*this));
	}

};
//end class VCluster

} /* end namespace xolotlCore */

#endif
