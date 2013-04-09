#ifndef REACTANT_H
#define REACTANT_H
//Begin section for file reactant.h
//TODO: Add definitions that you want preserved
//End section for file reactant.h




//<p>A Reactant is a general reacting body in a reaction network. It represents any body whose population can change with time due to reactions of any type.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class Reactant
{

    //Begin section for Reactant
    //TODO: Add attributes that you want preserved
    //End section for Reactant

    private:


        //<p>The total concentration of this Reactant.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double concentration;



        //<p>The name of this Reactant.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * name;



        //<p>An integer identification number for this reactant. </p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int id;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Reactant(Reactant & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~Reactant(); 



        //<p>The constructor.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Reactant(); 



        //<p>An alternative constructor that can be used to create a reactant with an initial concentration.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Reactant(double conc); 



        //<p>This operation returns the current concentration.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double getConcentration(); 



        //<p>This operation increases the concentration of the reactant by the specified amount.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void increaseConcentration(double deltaConc); 



        //<p>This operation decreases the concentration of the reactant by the specified amount.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void decreaseConcentration(double deltaConc); 



        //<p>This operation sets the concentration of the reactant to the specified amount.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void setConcentration(double conc); 



        //<p>This operation sets the concentration of the reactant to zero.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void zero(); 



        //<p>This operation returns the total flux of this reactant in the current network.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double getTotalFlux(); 



        //<p>This operation sets the collection of other reactants that make up the reaction network in which this reactant exists. </p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void setReactionNetwork(Reactant * * network); 



        //<p>This operation writes the contents of the reactant to a string.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * toString(); 



        //<p>This operation reads the contents of the reactant from a string.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void fromString(const char * input); 



        //<p>This operation writes the contents of the reactant to a string in XML.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * toXML(); 



        //<p>This operation reads the contents of the reactant from a string in XML.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void fromXML(const char * input); 



};  //end class Reactant



#endif
