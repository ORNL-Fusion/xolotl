#ifndef IDATAPROVIDER_H
#define IDATAPROVIDER_H
//Begin section for file idataprovider.h
//TODO: Add definitions that you want preserved
//End section for file idataprovider.h
#include "../vector.h"
#include "point.h"





//<p>IDataProvider describes the structure needed as a link between Xololt outputs and quantities that the user wants to plot. </p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class IDataProvider
{

    //Begin section for IDataProvider
    //TODO: Add attributes that you want preserved
    //End section for IDataProvider

    private:


        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Point * point;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual ~IDataProvider(); 



        //<p>Returns a string the name of the data.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getDataName() = 0; 



        //<p>Returns a string of the data unit.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getDataUnit() = 0; 



        //<p>Returns a collection of the data points.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual Vector getDataPoints() = 0; 



        //<p>Returns the value of the mean of all the data points.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual double getDataMean() = 0; 



        //<p>Method setting the name of the data.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setName(const char * name) = 0; 



        //<p>Method setting the unit of the data.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setUnits(const char * unit) = 0; 



        //<p>Method filling the data collection.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setPoints(Point Points) = 0; 



        //<p>Method returning the vector of quantity that will be plotted on the X axis. Quantity being x, y, z, t, or value.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual Vector getXVector() = 0; 



        //<p>Method returning the vector of quantity that will be plotted on the Y axis. Quantity being x, y, z, t, or value.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual Vector getYVector() = 0; 



        //<p>Method returning the vector of quantity that will be plotted on the Z axis. Quantity being x, y, z, t, or value.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual Vector getZVector() = 0; 



        //<p>Method returning the vector of time steps that will be used for the VideoPlots.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual Vector getTimeVector() = 0; 



};  //end class IDataProvider



#endif
