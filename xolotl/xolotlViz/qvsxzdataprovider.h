#ifndef QVSXZDATAPROVIDER_H
#define QVSXZDATAPROVIDER_H
//Begin section for file qvsxzdataprovider.h
//TODO: Add definitions that you want preserved
//End section for file qvsxzdataprovider.h
#include "dataprovider.h"
#include "../vector.h"





//<p>Subclass of DataProvider that will provide the methods to give the value, and X and Z data to a SurfacePlot.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class QvsXZDataProvider : DataProvider
{

    //Begin section for QvsXZDataProvider
    //TODO: Add attributes that you want preserved
    //End section for QvsXZDataProvider



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        QvsXZDataProvider(QvsXZDataProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~QvsXZDataProvider(); 



        //<p>Method returning a vector containing the 'x' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getXVector(); 



        //<p>Method returning a vector containing the 'z' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getYVector(); 



        //<p>Method returning a vector containing the 'value' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getZVector(); 



};  //end class QvsXZDataProvider



#endif
