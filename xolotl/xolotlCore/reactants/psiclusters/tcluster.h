#ifndef TCLUSTER_H
#define TCLUSTER_H

// Includes
#include "hcluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of tritium.
class TCluster : HCluster {

private:

	//! The number of tritium atoms in this cluster.
	int numT;

public:

	//! The constructor. All TClusters must be initialized with a size.
	TCluster(int nT);

	//! Destructor
	~TCluster();

};
//end class TCluster

} /* end namespace xolotlCore */
#endif
