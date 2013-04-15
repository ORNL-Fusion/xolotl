#ifndef VCLUSTER_H
#define VCLUSTER_H
#include "psicluster.h"

// This class represents a cluster composed entirely of atomic vacancies. 
class VCluster : PSICluster
{

    private:

        // The number of atomic vacancies in this cluster. 
        int numVac;

    public:

        VCluster(VCluster & arg); 

        ~VCluster(); 

        VCluster(); 

        ~VCluster(); 

        // The constructor. All VClusters must be initialized with a size. 
        VCluster(int nV); 

};  //end class VCluster

#endif
