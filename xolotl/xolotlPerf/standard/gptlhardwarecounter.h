#ifndef GPTLHARDWARECOUNTER_H
#define GPTLHARDWARECOUNTER_H
//Begin section for file gptlhardwarecounter.h
//TODO: Add definitions that you want preserved
//End section for file gptlhardwarecounter.h
#include "hardwarecounter.h"





//<p>The GPTLHardwareCounter class is instantiated by the StandardHandlerRegistry class and realizes the HardwareCounter interface to gather hardware performance counter data found by utilizing the PAPI (Performance Application Programming Interface) library via the General Purpose Timing Library (GPTL).</p><p></p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class GPTLHardwareCounter : HardwareCounter
{

    //Begin section for GPTLHardwareCounter
    //TODO: Add attributes that you want preserved
    //End section for GPTLHardwareCounter



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        GPTLHardwareCounter(GPTLHardwareCounter & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~GPTLHardwareCounter(); 



};  //end class GPTLHardwareCounter



#endif
