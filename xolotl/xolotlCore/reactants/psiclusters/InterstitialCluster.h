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

	/**
	 * This operation overrides getConnectivity() from the base class to
	 * provide the proper connectivity for an interstitial cluster.
	 * @return The connectivity array
	 */
    virtual std::vector<int> getConnectivity();

};
//end class InterstitialCluster

} /* end namespace xolotlCore */
#endif
