#ifndef QVSTIMEDATAPROVIDER_H
#define QVSTIMEDATAPROVIDER_H
//Begin section for file qvstimedataprovider.h
//TODO: Add definitions that you want preserved
//End section for file qvstimedataprovider.h
#include "dataprovider.h"
#include "../vector.h"





//<p>Subclass of DataProvider that will provide the methods to give the value and time data to a ScatterPlot.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class QvsTimeDataProvider : DataProvider
{

    //Begin section for QvsTimeDataProvider
    //TODO: Add attributes that you want preserved
    //End section for QvsTimeDataProvider



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        QvsTimeDataProvider(QvsTimeDataProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~QvsTimeDataProvider(); 



        //<p>Method returning a vector containing the 't' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getXVector(); 



        //<p>Method returning a vector containing the 'value' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getYVector(); 



};  //end class QvsTimeDataProvider



#endif
