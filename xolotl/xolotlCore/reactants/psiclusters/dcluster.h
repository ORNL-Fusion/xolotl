#ifndef DCLUSTER_H
#define DCLUSTER_H
//Begin section for file dcluster.h
//TODO: Add definitions that you want preserved
//End section for file dcluster.h
#include "hcluster.h"





//<p>This class represents a cluster composed entirely of deuterium.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class DCluster : HCluster
{

    //Begin section for DCluster
    //TODO: Add attributes that you want preserved
    //End section for DCluster

    private:


        //<p>The number of deuterium atoms in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numD;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        DCluster(DCluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~DCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~DCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        DCluster(); 



        //<p>The constructor. All DClusters must be initialized with a size.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        DCluster(int nH); 



};  //end class DCluster



#endif
