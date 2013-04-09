#ifndef TCLUSTER_H
#define TCLUSTER_H
//Begin section for file tcluster.h
//TODO: Add definitions that you want preserved
//End section for file tcluster.h
#include "hcluster.h"





//<p>This class represents a cluster composed entirely of tritium.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class TCluster : HCluster
{

    //Begin section for TCluster
    //TODO: Add attributes that you want preserved
    //End section for TCluster

    private:


        //<p>The number of tritium atoms in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numT;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        TCluster(TCluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~TCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~TCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        TCluster(); 



        //<p>The constructor. All TClusters must be initialized with a size.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        TCluster(int nH); 



};  //end class TCluster



#endif
