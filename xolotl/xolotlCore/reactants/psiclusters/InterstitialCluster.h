#ifndef INTERSTITIALCLUSTER_H
#define INTERSTITIALCLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of interstitial defects.
class InterstitialCluster : public PSICluster {

private:

public:

	/**
	 * The constructor. All InterstitialClusters must be initialized with
	 * a size.
	 */
	InterstitialCluster(int nI);

	//! Destructor
	~InterstitialCluster();

};
//end class InterstitialCluster

} /* end namespace xolotlCore */
#endif
