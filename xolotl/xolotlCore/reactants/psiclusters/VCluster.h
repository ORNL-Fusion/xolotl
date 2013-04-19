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

    public:

        /**
         * The constructor. All VClusters must be initialized with a size.
         * @param nV the number of atomic vacancies in the cluster
         */
        VCluster(int nV);

        //! Destructor
        ~VCluster(); 

};  //end class VCluster

} /* end namespace xolotlCore */

#endif
