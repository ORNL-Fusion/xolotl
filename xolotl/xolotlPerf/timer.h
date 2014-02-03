#ifndef TIMER_H
#define TIMER_H
//Begin section for file timer.h
//TODO: Add definitions that you want preserved
//End section for file timer.h




//<p>Realizations of this interface are responsible for the collection of performance timing statistics.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class Timer
{

    //Begin section for Timer
    //TODO: Add attributes that you want preserved
    //End section for Timer



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual ~Timer(); 



        //<p>Timer constructor that takes the argument name to distinguish specific timer.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Timer(const char * name); 



        //<p>This operations starts the Timer.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual void start() = 0; 



        //<p>This operation stops the Timer.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual void stop() = 0; 



        //<p>This operation returns the value of the Timer.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual long getValue() = 0; 



        //<p>This operation returns the units of the Timer.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual long getUnits() = 0; 



        //<p>This operation returns the name of the Timer.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getName() = 0; 



};  //end class Timer



#endif
