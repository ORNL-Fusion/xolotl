#ifndef QVSYZDATAPROVIDER_H
#define QVSYZDATAPROVIDER_H
//Begin section for file qvsyzdataprovider.h
//TODO: Add definitions that you want preserved
//End section for file qvsyzdataprovider.h
#include "dataprovider.h"
#include "../vector.h"





//<p>Subclass of DataProvider that will provide the methods to give the value, and Y and Z data to a SurfacePlot.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class QvsYZDataProvider : DataProvider
{

    //Begin section for QvsYZDataProvider
    //TODO: Add attributes that you want preserved
    //End section for QvsYZDataProvider



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        QvsYZDataProvider(QvsYZDataProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~QvsYZDataProvider(); 



        //<p>Method returning a vector containing the 'y' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getXVector(); 



        //<p>Method returning a vector containing the 'z' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getYVector(); 



        //<p>Method returning a vector containing the 'value' field of the collection of Point of the DataProvider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector getZVector(); 



};  //end class QvsYZDataProvider



#endif
