#ifndef QVSXZTIMEDATAPROVIDER_H
#define QVSXZTIMEDATAPROVIDER_H
//Begin section for file qvsxztimedataprovider.h
//TODO: Add definitions that you want preserved
//End section for file qvsxztimedataprovider.h
#include "dataprovider.h"
#include "../vector.h"





//<p>Subclass of DataProvider that will provide the methods to give the value, X and Z data, and time to a VideoPlot.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class QvsXZTimeDataProvider : DataProvider
{

    //Begin section for QvsXZTimeDataProvider
    //TODO: Add attributes that you want preserved
    //End section for QvsXZTimeDataProvider



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        QvsXZTimeDataProvider(QvsXZTimeDataProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~QvsXZTimeDataProvider(); 



        //<p>Method returning a vector containing the 'x' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getXVector(); 



        //<p>Method returning a vector containing the 'z' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getYVector(); 



        //<p>Method returning a vector containing the 'value' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getZVector(); 



        //<p>Method returning a vector containing the 'time' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getTimeVector(); 



};  //end class QvsXZTimeDataProvider



#endif
