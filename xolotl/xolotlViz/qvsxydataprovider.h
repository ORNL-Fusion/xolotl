#ifndef QVSXYDATAPROVIDER_H
#define QVSXYDATAPROVIDER_H
//Begin section for file qvsxydataprovider.h
//TODO: Add definitions that you want preserved
//End section for file qvsxydataprovider.h
#include "dataprovider.h"
#include "../vector.h"





//<p>Subclass of DataProvider that will provide the methods to give the value, and X and Y data to a SurfacePlot.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class QvsXYDataProvider : DataProvider
{

    //Begin section for QvsXYDataProvider
    //TODO: Add attributes that you want preserved
    //End section for QvsXYDataProvider



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        QvsXYDataProvider(QvsXYDataProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~QvsXYDataProvider(); 



        //<p>Method returning a vector containing the 'x' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getXVector(); 



        //<p>Method returning a vector containing the 'y' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getYVector(); 



        //<p>Method returning a vector containing the 'value' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getZVector(); 



};  //end class QvsXYDataProvider



#endif
