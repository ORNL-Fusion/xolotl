#ifndef DATAPROVIDERTESTER_H
#define DATAPROVIDERTESTER_H
//Begin section for file dataprovidertester.h
//TODO: Add definitions that you want preserved
//End section for file dataprovidertester.h
#include "../../visualization/dataprovider.h"





//<p>Class testing the data provider.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class DataProviderTester
{

    //Begin section for DataProviderTester
    //TODO: Add attributes that you want preserved
    //End section for DataProviderTester



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        DataProviderTester(DataProviderTester & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~DataProviderTester(); 



        //<p>Method checking you can set and get the unit of the data.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int checkUnits(); 



        //<p>Method checking you can add points to the data, and getting the data.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int checkPoints(); 



        //<p>Method checking you can set and get the name of the data.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int checkName(); 



        //<p>Method checking getting the mean of the value of the data points.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int checkGetDataMean(); 



};  //end class DataProviderTester



#endif
