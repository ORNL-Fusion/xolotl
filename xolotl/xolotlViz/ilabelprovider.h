#ifndef ILABELPROVIDER_H
#define ILABELPROVIDER_H
//Begin section for file ilabelprovider.h
//TODO: Add definitions that you want preserved
//End section for file ilabelprovider.h




//<p>ILabelProvider describes the structure needed to provide axis labels, titles, etc. for the plots.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class ILabelProvider
{

    //Begin section for ILabelProvider
    //TODO: Add attributes that you want preserved
    //End section for ILabelProvider



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual ~ILabelProvider(); 



        //<p>Returns the label of the X axis which doesn't necessary correspond to the physical X direction.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getXAxisLabel() = 0; 



        //<p>Returns the label of the Y axis which doesn't necessary correspond to the physical Y direction.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getYAxisLabel() = 0; 



        //<p>Returns the label of the Z axis which doesn't necessary correspond to the physical Z direction.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getZAxisLabel() = 0; 



        //<p>Returns a string of the time label.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getTimeAxisLabel() = 0; 



        //<p>Returns a string of the title label.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getTitleLabel() = 0; 



        //<p>Sets the name of the X axis of the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setXAxisLabel(const char * xLabel) = 0; 



        //<p>Sets the name of the Y axis of the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setYAxisLabel(const char * yLabel) = 0; 



        //<p>Sets the name of the Z axis of the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setZAxisLabel(const char * zLabel) = 0; 



        //<p>Sets the name of the time axis of the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setTimeAxisLabel(const char * timeLabel) = 0; 



        //<p>Sets the name of the title of the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setTitleLabel(const char * title) = 0; 



};  //end class ILabelProvider



#endif
