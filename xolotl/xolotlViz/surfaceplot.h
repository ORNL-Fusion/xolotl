#ifndef SURFACEPLOT_H
#define SURFACEPLOT_H
//Begin section for file surfaceplot.h
//TODO: Add definitions that you want preserved
//End section for file surfaceplot.h
#include "plot.h"
#include "../vector.h"
#include "point.h"





//<p>Plot the data value as a function of two different dimensions. The available PlottingStyle are POINTS, LINE, COLORMAP, and SURFACE.  It can be associated to QvsXYDataProvider, QvsXZDataProvider, or QvsYZDataProvider.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class SurfacePlot : Plot
{

    //Begin section for SurfacePlot
    //TODO: Add attributes that you want preserved
    //End section for SurfacePlot

    private:


        //<p>Vector containing three fields: one for the value and two for the directions as a function of which the value is plotted.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector data;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        SurfacePlot(SurfacePlot & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~SurfacePlot(); 



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



        //<p>Method getting the Z axis label with the help of the label provider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * getZAxisLabel(); 



};  //end class SurfacePlot



#endif
