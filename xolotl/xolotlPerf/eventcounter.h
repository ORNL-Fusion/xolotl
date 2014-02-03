#ifndef EVENTCOUNTER_H
#define EVENTCOUNTER_H
//Begin section for file eventcounter.h
//TODO: Add definitions that you want preserved
//End section for file eventcounter.h




//<p>Realizations of this interface are responsible for the collection of event performance counter data.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class EventCounter
{

    //Begin section for EventCounter
    //TODO: Add attributes that you want preserved
    //End section for EventCounter



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual ~EventCounter(); 



        //<p>Constructor that takes the name of the EventCounter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        EventCounter(const char * name); 



        //<p>This operation returns the value of the EventCounter, the frequency of the specified event.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual long getValue() = 0; 



        //<p>This operation returns the name of the EventCounter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getName() = 0; 



        //<p>This operation increments the EventCounter.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual void increment() = 0; 



};  //end class EventCounter



#endif
