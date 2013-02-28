#ifndef VCLUSTER_H
#define VCLUSTER_H
//Begin section for file vcluster.h
//TODO: Add definitions that you want preserved
//End section for file vcluster.h
#include "psicluster.h"





//<p>This class represents a cluster composed entirely of atomic vacancies.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class VCluster : PSICluster
{
    //Begin section for VCluster
    //TODO: Add attributes that you want preserved
    //End section for VCluster

    private:


        //<p>The number of atomic vacancies in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numVac;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        VCluster(VCluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~VCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        VCluster(VCluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~VCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        VCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~VCluster(); 



        //get numVac
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int & get_numVac(); 



        //set numVac
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void set_numVac(int & numVac); 



        //<p>The constructor. All VClusters must be initialized with a size.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        VCluster(int nV); 



};  //end class VCluster



#endif
