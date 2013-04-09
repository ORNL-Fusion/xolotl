#ifndef INTERSTITIALCLUSTER_H
#define INTERSTITIALCLUSTER_H
//Begin section for file interstitialcluster.h
//TODO: Add definitions that you want preserved
//End section for file interstitialcluster.h
#include "psicluster.h"





//<p>This class represents a cluster composed entirely of interstitial defects.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class InterstitialCluster : PSICluster
{

    //Begin section for InterstitialCluster
    //TODO: Add attributes that you want preserved
    //End section for InterstitialCluster

    private:


        //<p>The number of interstitial defects in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numInt;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        InterstitialCluster(InterstitialCluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~InterstitialCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        InterstitialCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~InterstitialCluster(); 



        //<p>The constructor. All InterstitialClusters must be initialized with a size.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        InterstitialCluster(int nI); 



};  //end class InterstitialCluster



#endif
