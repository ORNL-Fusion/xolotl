#ifndef HECLUSTER_H
#define HECLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of helium.
 */
class HeCluster: public PSICluster {

private:

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	HeCluster():PSICluster(1) {}

public:

	/**
	 * The constructor. All HeClusters must be initialized with a size.
	 * @param nHe the number of helium atoms in the cluster
	 */
	HeCluster(int nHe);

	/**
	 * Destructor
	 */
	~HeCluster();

	/**
	 * This operation overrides getConnectivity() from the base class to
	 * provide the proper connectivity for a He cluster.
	 * @return The connectivity array
	 */
    virtual std::vector<int> getConnectivity();

};
//end class HeCluster

} /* namespace xolotlCore */
#endif
