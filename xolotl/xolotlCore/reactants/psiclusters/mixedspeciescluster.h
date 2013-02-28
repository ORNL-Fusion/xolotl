#ifndef MIXEDSPECIESCLUSTER_H
#define MIXEDSPECIESCLUSTER_H
//Begin section for file mixedspeciescluster.h
//TODO: Add definitions that you want preserved
//End section for file mixedspeciescluster.h
#include "psicluster.h"





//<p>This class represents a cluster composed of multiple atomic species, vacancies and interstitials.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class MixedSpeciesCluster : PSICluster
{
    //Begin section for MixedSpeciesCluster
    //TODO: Add attributes that you want preserved
    //End section for MixedSpeciesCluster

    private:


        //<p>The total number of atoms, vacancies and interstitials making up this cluster across all species.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int size;



        //<p>The number of hydrogen atoms in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numH;



        //<p>The number of interstitial defects in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numInt;



        //<p>The number of deuterium atoms in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numD;



        //<p>The number of helium atoms in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numHe;



        //<p>The number of tritium atoms in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numT;



        //<p>The number of atomic vacancies in this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int numVac;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        MixedSpeciesCluster(MixedSpeciesCluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~MixedSpeciesCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        MixedSpeciesCluster(MixedSpeciesCluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~MixedSpeciesCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        MixedSpeciesCluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~MixedSpeciesCluster(); 



        //get size
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int & get_size(); 



        //set size
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void set_size(int & size); 



        //get numH
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int & get_numH(); 



        //set numH
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void set_numH(int & numH); 



        //get numInt
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int & get_numInt(); 



        //set numInt
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void set_numInt(int & numInt); 



        //get numD
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int & get_numD(); 



        //set numD
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void set_numD(int & numD); 



        //get numHe
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int & get_numHe(); 



        //set numHe
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void set_numHe(int & numHe); 



        //get numT
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int & get_numT(); 



        //set numT
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void set_numT(int & numT); 



        //get numVac
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int & get_numVac(); 



        //set numVac
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void set_numVac(int & numVac); 



        //<p>The constructor. All MixedSpeciesClusters must be initialized with two arrays that describe the species of which it is composed. The arrays contain the names of the species and the respective concentration of each species. The species must be one of {H,He,I,V,D,T}.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        MixedSpeciesCluster(int numPerSpecies, const char * speciesNames); 



        //<p>This operation returns the total generation rate due to emission for this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double getGenByEm(); 



        //<p>This operation returns the total generation rate due to annihilation for this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double getAnnByEm(); 



        //<p>This operation returns the number of a given "species" within this cluster by passing one of {H,He,I,V,D,T} as an input argument.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int getSpeciesSize(const char * speciesName); 



};  //end class MixedSpeciesCluster



#endif
