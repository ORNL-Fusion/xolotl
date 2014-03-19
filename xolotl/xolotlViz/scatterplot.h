#ifndef SCATTERPLOT_H
#define SCATTERPLOT_H
//Begin section for file scatterplot.h
//TODO: Add definitions that you want preserved
//End section for file scatterplot.h
#include "plot.h"
#include "../vector.h"
#include "point.h"





//<p>Plot the data value as a function of one dimension. Available PlottingStyle are POINTS or LINE. It can be associated to QvsXDataProvider, QvsYDataProvider, QvsZDataProvider, or QvsTimeDataProvider.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class ScatterPlot : Plot
{

    //Begin section for ScatterPlot
    //TODO: Add attributes that you want preserved
    //End section for ScatterPlot

    private:


        //<p>Vector containing two fields: one for the value and one for the direction as a function of which the value is plotted.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector data;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ScatterPlot(ScatterPlot & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~ScatterPlot(); 



        //<p>Method managing everything that is related to the rendering of a plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int render(bool isCumulative=false); 



        //<p>Method returning the data points that are stored in the data vector.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Point getPoints(); 



        //<p>Method getting the X axis label with the help of the label provider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * getXAxisLabel(); 



        //<p>Method getting the Y axis label with the help of the label provider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * getYAxisLabel(); 



};  //end class ScatterPlot



#endif
