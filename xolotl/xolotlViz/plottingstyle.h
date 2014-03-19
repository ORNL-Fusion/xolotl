#ifndef PLOTTINGSTYLE_H
#define PLOTTINGSTYLE_H
//Begin section for file plottingstyle.h
//TODO: Add definitions that you want preserved
//End section for file plottingstyle.h


//<p>Enumeration of the different plotting styles available.</p>

//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
enum PlottingStyle
{



        //<p>The data will be shown as points, whether it is a ScatterPlot, a SurfacePlot, or a VideoPlot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        POINTS,


        //<p>Lines between the data points will be drawn for a ScatterPlot. A grid between the data points will be represented for a SurfacePlot or a VideoPlot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        LINE,


        //<p>Only available for SurfacePlot and VideoPlot. The value at each data point will be represented by a color.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        COLORMAP,


        //<p>Only available for SurfacePlot and VideoPlot. The rendering will be a 3D surface as well as colors.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        SURFACE

}; //end enum PlottingStyle 


#endif
