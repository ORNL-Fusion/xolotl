#ifndef QVSXDATAPROVIDER_H
#define QVSXDATAPROVIDER_H
//Begin section for file qvsxdataprovider.h
//TODO: Add definitions that you want preserved
//End section for file qvsxdataprovider.h
#include "dataprovider.h"
#include "../vector.h"





//<p>Subclass of DataProvider that will provide the methods to give the value and X data to a ScatterPlot.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class QvsXDataProvider : DataProvider
{

    //Begin section for QvsXDataProvider
    //TODO: Add attributes that you want preserved
    //End section for QvsXDataProvider



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        QvsXDataProvider(QvsXDataProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~QvsXDataProvider(); 



        //<p>Method returning a vector containing the 'x' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getXVector(); 



        //<p>Method returning a vector containing the 'value' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getYVector(); 



};  //end class QvsXDataProvider



#endif
