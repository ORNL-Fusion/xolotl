#ifndef DUMMYHANDLERREGISTRY_H
#define DUMMYHANDLERREGISTRY_H
//Begin section for file dummyhandlerregistry.h
//TODO: Add definitions that you want preserved
//End section for file dummyhandlerregistry.h
#include "ihandlerregistry.h"
#include "dummytimer.h" //Dependency Generated Source:DummyHandlerRegistry Target:DummyTimer
#include "dummyeventcounter.h" //Dependency Generated Source:DummyHandlerRegistry Target:DummyEventCounter
#include "dummyhardwarecounter.h" //Dependency Generated Source:DummyHandlerRegistry Target:DummyHardwareCounter





//<p>The DummyHandlerRegistry class realizes the interface IHandlerRegistry to support the "always on" performance data collection acquired by implementing the performance interfaces Timer, EventCounter, and HardwareCounter.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class DummyHandlerRegistry : IHandlerRegistry
{

    //Begin section for DummyHandlerRegistry
    //TODO: Add attributes that you want preserved
    //End section for DummyHandlerRegistry



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        DummyHandlerRegistry(DummyHandlerRegistry & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~DummyHandlerRegistry(); 



};  //end class DummyHandlerRegistry



#endif
