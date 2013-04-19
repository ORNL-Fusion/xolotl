#ifndef HCLUSTER_H
#define HCLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of hydrogen.
class HCluster : public PSICluster
{

    private:

    public:

    	/**
    	 * The constructor. All HClusters must be initialized with a size.
    	 * @param nH the number of hydrogen atoms in the cluster
    	 */
        HCluster(int nH);

        //! The destructor
        ~HCluster(); 

};  //end class HCluster

} /* end namespace xolotlCore */
#endif
