#ifndef DATAPROVIDER_H
#define DATAPROVIDER_H
//Begin section for file dataprovider.h
//TODO: Add definitions that you want preserved
//End section for file dataprovider.h
#include "idataprovider.h"
#include "../vector.h"





//<p>Realization of the IDataProvider interface. This is a general class with general methods, to actually get data from the data provider, one needs to use the subclasses.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class DataProvider : IDataProvider
{

    //Begin section for DataProvider
    //TODO: Add attributes that you want preserved
    //End section for DataProvider

    private:


        //<p>Collection of data points.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector data;



        //<p>Name of the data contained in the data attribute.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * dataName;



        //<p>Unit of the data contained in the data attribute.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * dataUnit;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        DataProvider(DataProvider & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~DataProvider(); 



        //<p>Virtual method returning the vector of quantity that will be plotted on the X axis. Quantity being x, y, z, t, or value.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual Vector getXVector() = 0; 



        //<p>Virtual method returning the vector of quantity that will be plotted on the Y axis. Quantity being x, y, z, t, or value.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual Vector getYVector() = 0; 



        //<p>Virtual method returning the vector of quantity that will be plotted on the Z axis. Quantity being x, y, z, t, or value.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual Vector getZVector() = 0; 



        //<p>Virtual method returning the vector of time steps that will be used for the VideoPlots.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual Vector getTimeVector() = 0; 



};  //end class DataProvider



#endif
