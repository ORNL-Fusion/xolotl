#ifndef INTERSTITIALCLUSTER_H
#define INTERSTITIALCLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of interstitial defects.
class InterstitialCluster: public PSICluster {

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

	/**
	 * This operation returns the total change in this cluster due to
	 * dissociation.
	 * @param temperature The temperature at which to calculate the Diffusion Coefficient
	 * @return The flux due to dissociation.
	 */
	virtual double getDissociationFlux(const double temperature);

	/**
	 * This virtual method is for subclasses to specialize
	 * to return their representative cluster map, which is a mapping
	 * of which species exist in the cluster to the integer number
	 * of each species.
	 *
	 * @return
	 */
	virtual std::map<std::string, int> getClusterMap();

};
//end class InterstitialCluster

} /* end namespace xolotlCore */
#endif
