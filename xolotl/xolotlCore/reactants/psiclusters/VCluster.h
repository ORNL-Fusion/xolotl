#ifndef VCLUSTER_H
#define VCLUSTER_H

// Includes
#include "PSICluster.h"

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
			PSICluster(1) {
	}

public:

	/**
	 * The constructor. All VClusters must be initialized with a size.
	 * @param nV the number of atomic vacancies in the cluster
	 */
	VCluster(int nV);

	//! Destructor
	~VCluster();

	/**
	 * This operation overrides getConnectivity() from the base class to
	 * provide the proper connectivity for a V cluster.
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
//end class VCluster

} /* end namespace xolotlCore */

#endif
