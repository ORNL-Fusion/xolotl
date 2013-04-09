#ifndef HCLUSTER_H
#define HCLUSTER_H
//Begin section for file hcluster.h
//TODO: Add definitions that you want preserved
//End section for file hcluster.h
#include "psicluster.h"





//<p>This class represents a cluster composed entirely of hydrogen.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class HCluster : PSICluster
{

    //Begin section for HCluster
    //TODO: Add attributes that you want preserved
    //End section for HCluster

    private:


        //<p>The number of hydrogen atoms in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numH;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        HCluster(HCluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~HCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~HCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        HCluster(); 



        //<p>The constructor. All HClusters must be initialized with a size.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        HCluster(int nH); 



};  //end class HCluster



#endif
