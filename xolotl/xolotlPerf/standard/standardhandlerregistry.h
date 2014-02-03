#ifndef STANDARDHANDLERREGISTRY_H
#define STANDARDHANDLERREGISTRY_H
//Begin section for file standardhandlerregistry.h
//TODO: Add definitions that you want preserved
//End section for file standardhandlerregistry.h
#include "ihandlerregistry.h"
#include "gptltimer.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLTimer
#include "gptleventcounter.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLEventCounter
#include "gptlhardwarecounter.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLHardwareCounter





//<p>The StandardHandlerRegistry class realizes the interface IHandlerRegistry to acquire the performance data found by implementing the performance interfaces Timer, EventCounter, and HardwareCounter.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class StandardHandlerRegistry : IHandlerRegistry
{

    //Begin section for StandardHandlerRegistry
    //TODO: Add attributes that you want preserved
    //End section for StandardHandlerRegistry



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        StandardHandlerRegistry(StandardHandlerRegistry & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~StandardHandlerRegistry(); 



};  //end class StandardHandlerRegistry



#endif
