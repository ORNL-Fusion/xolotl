#ifndef PSICLUSTER_H
#define PSICLUSTER_H
//Begin section for file psicluster.h
//TODO: Add definitions that you want preserved
//End section for file psicluster.h
#include "reactant.h"





//<p>The PSICluster class is a Reactant that is specialized to work for simulations of plasma-surface interactions. It provides special routines for calculating the total flux due to production and dissociation and obtaining the cluster size.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class PSICluster : Reactant
{

    //Begin section for PSICluster
    //TODO: Add attributes that you want preserved
    //End section for PSICluster

    private:


        //<p>The diffusion factor, D_0, that is used to calculate the diffusion coefficient for this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double diffusionFactor;



        //<p>The migration energy for this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double migrationEnergy;



        //<p>The binding energy for this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double bindingEnergy;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        PSICluster(PSICluster & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~PSICluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        PSICluster(); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~PSICluster(); 



        //<p>This operation returns the total change in this cluster due to dissociation.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual double getDissociationFlux() = 0; 



        //<p>This operation returns the total change in this cluster due to production.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual double getProductionFlux() = 0; 



        //<p>This operation returns the total size of the cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int getSize() = 0; 



        //<p>This operation returns the total generation rate of this cluster due to captures.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual double getGenByCapt() = 0; 



        //<p>This operation returns the total generation rate of this cluster due to annihilation.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual double getGenByAnn() = 0; 



        //<p>This operation retrieves the diffusion factor, D_0, that is used to calculate the diffusion coefficient for this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double getDiffusionFactor(); 



        //<p>This operation retrieves the migration energy for this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double getMigrationEnergy(); 



        //<p>This operation retrieves the binding energy for this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double getBindingEnergy(); 



        //<p>This operation returns the diffusion coefficient for this cluster and is based on the diffusion factor.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual double getDiffusionCoefficient() = 0; 



        //<p>This operation sets the diffusion factor, D_0, that is used to calculate the diffusion coefficient for this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void setDiffusionFactor(double factor); 



        //<p>This operation sets the migration energy for this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void setMigrationEnergy(double em); 



        //<p>This operation sets the binding energy for this cluster.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void setBindingEnergy(double eb); 



};  //end class PSICluster



#endif
