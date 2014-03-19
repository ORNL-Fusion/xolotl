#ifndef QVSZDATAPROVIDER_H
#define QVSZDATAPROVIDER_H
//Begin section for file qvszdataprovider.h
//TODO: Add definitions that you want preserved
//End section for file qvszdataprovider.h
#include "dataprovider.h"
#include "../vector.h"





//<p>Subclass of DataProvider that will provide the methods to give the value and Z data to a ScatterPlot.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class QvsZDataProvider : DataProvider
{

    //Begin section for QvsZDataProvider
    //TODO: Add attributes that you want preserved
    //End section for QvsZDataProvider



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        QvsZDataProvider(QvsZDataProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~QvsZDataProvider(); 



        //<p>Method returning a vector containing the 'z' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getXVector(); 



        //<p>Method returning a vector containing the 'value' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getYVector(); 



};  //end class QvsZDataProvider



#endif
