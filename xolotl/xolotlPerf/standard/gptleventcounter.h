#ifndef GPTLEVENTCOUNTER_H
#define GPTLEVENTCOUNTER_H
//Begin section for file gptleventcounter.h
//TODO: Add definitions that you want preserved
//End section for file gptleventcounter.h
#include "eventcounter.h"





//<p>The GPTLEventCounter class is instantiated by the StandardHandlerRegistry class and realizes the EventCounter interface to access event performance counter data found via the General Purpose Timing Library (GPTL).</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class GPTLEventCounter : EventCounter
{

    //Begin section for GPTLEventCounter
    //TODO: Add attributes that you want preserved
    //End section for GPTLEventCounter



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        GPTLEventCounter(GPTLEventCounter & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~GPTLEventCounter(); 



};  //end class GPTLEventCounter



#endif
