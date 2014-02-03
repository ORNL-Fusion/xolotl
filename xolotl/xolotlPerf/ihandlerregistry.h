#ifndef IHANDLERREGISTRY_H
#define IHANDLERREGISTRY_H
//Begin section for file ihandlerregistry.h
//TODO: Add definitions that you want preserved
//End section for file ihandlerregistry.h
#include "timer.h" //Dependency Generated Source:IHandlerRegistry Target:Timer
#include "eventcounter.h" //Dependency Generated Source:IHandlerRegistry Target:EventCounter
#include "hardwarecounter.h" //Dependency Generated Source:IHandlerRegistry Target:HardwareCounter
#include "hardwarequantities.h"





//<p>Realizations of this interface are responsible for the collection of performance data.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class IHandlerRegistry
{

    //Begin section for IHandlerRegistry
    //TODO: Add attributes that you want preserved
    //End section for IHandlerRegistry



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual ~IHandlerRegistry(); 



        //<p>This operation returns the Timer specified by the parameter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual Timer getTimer(const char * name) = 0; 



        //<p>This operation returns the EventCounter specified by the parameter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual EventCounter getEventCounter(const char * name) = 0; 



        //<p>This operation returns the specified HardwareCounter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual HardwareCounter getHardwareCounter(const char * name,  * quantities) = 0; 



        //<p>This operation returns a list of values of the, initially specified, PAPI preset quantities monitored by the HardwareCounter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual HardwareQuantities * getHardwareQuantities() = 0; 



        //<p>This operation outputs the information gathered.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual void dump( out) = 0; 



};  //end class IHandlerRegistry



#endif
