#ifndef VIDEOPLOT_H
#define VIDEOPLOT_H
//Begin section for file videoplot.h
//TODO: Add definitions that you want preserved
//End section for file videoplot.h
#include "plot.h"
#include "../java/util/vector_0.h"
#include "point.h"



class SurfacePlot; //Dependency Generated Source:VideoPlot Target:SurfacePlot





//<p>Plot the data value as a function of two different spatial dimensions for each time step and change the time step with time to have a video-like rendering, each frame being a SurfacePlot. The available PlottingStyle are POINTS, LINE, COLORMAP, and SURFACE. It can be associated to QvsXYTimeDataProvider, QvsXZTimeDataProvider, or QvsYZTimeDataProvider.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class VideoPlot : Plot
{

    //Begin section for VideoPlot
    //TODO: Add attributes that you want preserved
    //End section for VideoPlot

    private:


        //<p>Vector containing four fields: one for the value, two for the spatial directions as a function of which the value is plotted, and the last one for the time step.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Vector_0 data;



        //<p>Number of frames shown per second.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double frameRate;



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        SurfacePlot * surfacePlot;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        VideoPlot(VideoPlot & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~VideoPlot(); 



        //<p>Method managing everything that is related to the rendering of a plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int render(bool isCumulative=false, double waitingTime=1); 



        //<p>Sets the frame rate for VideoPlot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int setFrameRate(double frameRate); 



        //<p>Get the frame rate from VideoPlot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double getFrameRate(); 



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



        //<p>Method getting the time axis label with the help of the label provider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * getTimeAxisLabel(); 



};  //end class VideoPlot



#endif
