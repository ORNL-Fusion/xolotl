#ifndef QVSYDATAPROVIDER_H
#define QVSYDATAPROVIDER_H
//Begin section for file qvsydataprovider.h
//TODO: Add definitions that you want preserved
//End section for file qvsydataprovider.h
#include "dataprovider.h"
#include "../vector.h"





//<p>Subclass of DataProvider that will provide the methods to give the value and Y data to a ScatterPlot.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class QvsYDataProvider : DataProvider
{

    //Begin section for QvsYDataProvider
    //TODO: Add attributes that you want preserved
    //End section for QvsYDataProvider



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        QvsYDataProvider(QvsYDataProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~QvsYDataProvider(); 



        //<p>Method returning a vector containing the 'y' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getXVector(); 



        //<p>Method returning a vector containing the 'value' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getYVector(); 



};  //end class QvsYDataProvider



#endif
