#ifndef STANDARDHANDLERREGISTRYTESTER_H
#define STANDARDHANDLERREGISTRYTESTER_H
//Begin section for file standardhandlerregistrytester.h
//TODO: Add definitions that you want preserved
//End section for file standardhandlerregistrytester.h
#include "timertester.h"
#include "hardwarecountertester.h"
#include "../../performance/standard/standardhandlerregistry.h"
#include "../../performance/timer.h"
#include "../../performance/eventcounter.h"
#include "../../performance/hardwarecounter.h"





//<p>This class is responsible for testing StandardHandlerRegistry.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class StandardHandlerRegistryTester
{

    //Begin section for StandardHandlerRegistryTester
    //TODO: Add attributes that you want preserved
    //End section for StandardHandlerRegistryTester

    private:


        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        TimerTester testTimer;



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
         testEventCounter;



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        HardwareCounterTester testHardwareCounter;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        StandardHandlerRegistryTester(StandardHandlerRegistryTester & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~StandardHandlerRegistryTester(); 



        //<p>This operation checks the Timer.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void checkTimer(); 



        //<p>This operation checks the EventCounter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void checkEventCounter(); 



        //<p>This operation checks the HardwareCounter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        void checkHardwareCounter(); 



};  //end class StandardHandlerRegistryTester



#endif
