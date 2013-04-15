#ifndef HECLUSTER_H
#define HECLUSTER_H

// Includes
#include "psicluster.h"

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of helium.
 */
class HeCluster: PSICluster {

private:

	//! The number of helium atoms in this cluster.
	int numHe;

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

};
//end class HeCluster

} /* namespace xolotlCore */
#endif
