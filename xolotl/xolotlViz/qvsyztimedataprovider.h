#ifndef QVSYZTIMEDATAPROVIDER_H
#define QVSYZTIMEDATAPROVIDER_H
//Begin section for file qvsyztimedataprovider.h
//TODO: Add definitions that you want preserved
//End section for file qvsyztimedataprovider.h
#include "dataprovider.h"
#include "../vector.h"





//<p>Subclass of DataProvider that will provide the methods to give the value, Y and Z data, and time to a VideoPlot.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class QvsYZTimeDataProvider : DataProvider
{

    //Begin section for QvsYZTimeDataProvider
    //TODO: Add attributes that you want preserved
    //End section for QvsYZTimeDataProvider



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        QvsYZTimeDataProvider(QvsYZTimeDataProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~QvsYZTimeDataProvider(); 



        //<p>Method returning a vector containing the 'y' field of the collection of Point of the DataProvider.</p>
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



};  //end class QvsYZTimeDataProvider



#endif
