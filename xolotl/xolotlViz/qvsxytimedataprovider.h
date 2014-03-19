#ifndef QVSXYTIMEDATAPROVIDER_H
#define QVSXYTIMEDATAPROVIDER_H
//Begin section for file qvsxytimedataprovider.h
//TODO: Add definitions that you want preserved
//End section for file qvsxytimedataprovider.h
#include "dataprovider.h"
#include "../vector.h"





//<p>Subclass of DataProvider that will provide the methods to give the value, X and Y data, and time to a VideoPlot.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class QvsXYTimeDataProvider : DataProvider
{

    //Begin section for QvsXYTimeDataProvider
    //TODO: Add attributes that you want preserved
    //End section for QvsXYTimeDataProvider



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        QvsXYTimeDataProvider(QvsXYTimeDataProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~QvsXYTimeDataProvider(); 



        //<p>Method returning a vector containing the 'x' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getXVector(); 



        //<p>Method returning a vector containing the 'y' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getYVector(); 



        //<p>Method returning a vector containing the 'value' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getZVector(); 



        //<p>Method returning a vector containing the 'time' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getTimeVector(); 



};  //end class QvsXYTimeDataProvider



#endif
