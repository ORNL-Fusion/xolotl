#ifndef LABELPROVIDER_H
#define LABELPROVIDER_H
//Begin section for file labelprovider.h
//TODO: Add definitions that you want preserved
//End section for file labelprovider.h
#include "ilabelprovider.h"





//<p>Class realizing the interface ILavelProvider. Contains different string attributes representing the labels for the plots. LabelProvider is a class attached to a single Plot at its creation.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class LabelProvider : ILabelProvider
{

    //Begin section for LabelProvider
    //TODO: Add attributes that you want preserved
    //End section for LabelProvider

    private:


        //<p>The label of the X axis of the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * xAxisLabel;



        //<p>The label of the Y axis of the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * yAxisLabel;



        //<p>The label of the Z axis of the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * zAxisLabel;



        //<p>The label for the time steps.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * timeAxisLabel;



        //<p>Title label for the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * titleLabel;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        LabelProvider(LabelProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~LabelProvider(); 



};  //end class LabelProvider



#endif
