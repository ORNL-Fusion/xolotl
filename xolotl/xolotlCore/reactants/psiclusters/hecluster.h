#ifndef HECLUSTER_H
#define HECLUSTER_H
//Begin section for file hecluster.h
//TODO: Add definitions that you want preserved
//End section for file hecluster.h
#include "psicluster.h"





//<p>This class represents a cluster composed entirely of helium.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class HeCluster : PSICluster
{
    //Begin section for HeCluster
    //TODO: Add attributes that you want preserved
    //End section for HeCluster

    private:


        //<p>The number of helium atoms in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numHe;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        HeCluster(HeCluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~HeCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        HeCluster(HeCluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~HeCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        HeCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~HeCluster(); 



        //get numHe
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int & get_numHe(); 



        //set numHe
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void set_numHe(int & numHe); 



        //<p>The constructor. All HeClusters must be initialized with a size.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        HeCluster(int nHe); 



};  //end class HeCluster



#endif
