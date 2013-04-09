#include "reactant.h"
//Begin section for file reactant.cpp
//TODO: Add definitions that you want preserved
//End section for file reactant.cpp


//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
Reactant::Reactant(Reactant & arg) : 
name("")
{
    //TODO Auto-generated method stub
    this->concentration = arg.concentration;
    this->id = arg.id;
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
Reactant::~Reactant() 
{
    //TODO Auto-generated method stub
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
Reactant::Reactant() : 
concentration(0.0),name(""),id(0)
{
    //TODO Auto-generated method stub
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
Reactant::Reactant(double conc) : 
concentration(0.0),name(""),id(0)
{
    //TODO Auto-generated method stub
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
double Reactant::getConcentration() 
{
    //TODO Auto-generated method stub
    return 0;
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
void Reactant::increaseConcentration(double deltaConc) 
{
    //TODO Auto-generated method stub
    return;
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
void Reactant::decreaseConcentration(double deltaConc) 
{
    //TODO Auto-generated method stub
    return;
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
void Reactant::setConcentration(double conc) 
{
    //TODO Auto-generated method stub
    return;
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
void Reactant::zero() 
{
    //TODO Auto-generated method stub
    return;
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
double Reactant::getTotalFlux() 
{
    //TODO Auto-generated method stub
    return 0;
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
void Reactant::setReactionNetwork(Reactant * * network) 
{
    //TODO Auto-generated method stub
    return;
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
const char * Reactant::toString() 
{
    //TODO Auto-generated method stub
    return 0;
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
void Reactant::fromString(const char * input) 
{
    //TODO Auto-generated method stub
    return;
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
const char * Reactant::toXML() 
{
    //TODO Auto-generated method stub
    return 0;
}
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
void Reactant::fromXML(const char * input) 
{
    //TODO Auto-generated method stub
    return;
}
