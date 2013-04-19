#ifndef DCLUSTER_H
#define DCLUSTER_H

// Includes
#include "HCluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of deuterium.
class DCluster : public HCluster {

private:

public:

	//! The constructor. All DClusters must be initialized with a size.
	DCluster(int nD);

	//! Destructor
	~DCluster();

};
//end class DCluster

} /* end namespace xolotlCore */
#endif
