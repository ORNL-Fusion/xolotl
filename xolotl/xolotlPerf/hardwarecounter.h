#ifndef HARDWARECOUNTER_H
#define HARDWARECOUNTER_H
//Begin section for file hardwarecounter.h
//TODO: Add definitions that you want preserved
//End section for file hardwarecounter.h
#include "hardwarequantities_0.h"





//<p>Realizations of this interface are responsible for the collection of hardware performance counter data.</p><p></p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class HardwareCounter
{

    //Begin section for HardwareCounter
    //TODO: Add attributes that you want preserved
    //End section for HardwareCounter



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual ~HardwareCounter(); 



        //<p>Constructor that takes the name of the HardwareCounter and a list of the different quantities that specific HardwareCounter should monitor.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        HardwareCounter(const char * name, HardwareQuantities_0 * quantities_); 



        //<p>This operation returns a list of values of the, initially specified, PAPI preset quantities monitored by the HardwareCounter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual long * getValues() = 0; 



        //<p>This operation returns the name of the HardwareCounter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getName() = 0; 



        //<p>This operation increments the HardwareCounter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual void increment() = 0; 



};  //end class HardwareCounter



#endif
