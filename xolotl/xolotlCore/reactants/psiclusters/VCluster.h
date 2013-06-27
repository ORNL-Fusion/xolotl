#ifndef VCLUSTER_H
#define VCLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of atomic vacancies.
 */
class VCluster : public PSICluster
{

    private:

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	VCluster():PSICluster(1) {}

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

};  //end class VCluster

} /* end namespace xolotlCore */

#endif
