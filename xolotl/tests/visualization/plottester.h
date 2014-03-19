#ifndef PLOTTESTER_H
#define PLOTTESTER_H
//Begin section for file plottester.h
//TODO: Add definitions that you want preserved
//End section for file plottester.h
#include "../../visualization/plot.h"





//<p>Class testing the general Plot class.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class PlotTester
{

    //Begin section for PlotTester
    //TODO: Add attributes that you want preserved
    //End section for PlotTester



    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        PlotTester(PlotTester & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~PlotTester(); 



        //<p>Method checking getting and setting the unit.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int checkUnit(); 



        //<p>Method checking the ability to choose a PlottingStyle.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int checkPlottingStyle(); 



        //<p>Method checking everything related the label provider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int checkLabelProvider(); 



        //<p>Method checking the ability to change the title.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int checkTitle(); 



        //<p>Method checking everything that is related to the legend.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int checkLegend(); 



        //<p>Method checking the writing of the file.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int checkWrite(); 



        //<p>Method checking everything related to the data provider in the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int checkDataProvider(); 



};  //end class PlotTester



#endif
