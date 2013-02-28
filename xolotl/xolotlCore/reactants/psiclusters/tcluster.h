#ifndef TCLUSTER_H
#define TCLUSTER_H
//Begin section for file tcluster.h
//TODO: Add definitions that you want preserved
//End section for file tcluster.h
#include "hcluster.h"
#include "psicluster.h"





//<p>This class represents a cluster composed entirely of tritium.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class TCluster : HCluster
{
    //Begin section for TCluster
    //TODO: Add attributes that you want preserved
    //End section for TCluster

    public:


        //<p>This class represents a cluster composed entirely of hydrogen.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        class HCluster : PSICluster
        {
            //Begin section for TCluster::HCluster
            //TODO: Add attributes that you want preserved
            //End section for TCluster::HCluster

            private:


                //<p>The number of hydrogen atoms in this cluster.</p>
                //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
                int numH;




            public:

                //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
                HCluster(TCluster::HCluster & arg); 



                //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
                ~HCluster(); 



                //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
                HCluster(TCluster::HCluster & arg); 



                //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
                ~HCluster(); 



                //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
                HCluster(); 



                //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
                ~HCluster(); 



                //get numH
                //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
                int & get_numH(); 



                //set numH
                //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
                void set_numH(int & numH); 



                //<p>The constructor. All HClusters must be initialized with a size.</p>
                //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
                HCluster(int nH); 



        };  //end class HCluster


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
        TCluster(TCluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~TCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        TCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~TCluster(); 



        //get numT
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int & get_numT(); 



        //set numT
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void set_numT(int & numT); 



        //<p>The constructor. All TClusters must be initialized with a size.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        TCluster(int nH); 



};  //end class TCluster



#endif
