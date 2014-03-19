#ifndef PLOT_H
#define PLOT_H
//Begin section for file plot.h
//TODO: Add definitions that you want preserved
//End section for file plot.h
#include "iplot.h"
#include "plottingstyle.h" //Dependency Generated Source:Plot Target:PlottingStyle
#include "dataprovider.h"
#include "labelprovider.h"





//<p>Plot is the class that realizes the interface Iplot. It is a general class that provides general methods, but to actual plot anything, the user needs to use one of its subclasses.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class Plot : IPlot
{

    //Begin section for Plot
    //TODO: Add attributes that you want preserved
    //End section for Plot

    private:


        //<p>Unit of the data that is being plotted.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * unit;



        //<p>Choice of PlottingStyle.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        PlottingStyle style;



        //<p>Label that will appear of the X axis of the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * XAxisLabel;



        //<p>Label that will appear of the Y axis of the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * YAxisLabel;



        //<p>Label that will appear of the Z axis of the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * ZAxisLabel;



        //<p>Title of the plot. Can be set by the user by using setTitle().</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * title;



        //<p>If it is equal to True, the legend will be displayed.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        bool enableLegend;



        //<p>Data provider used for the plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        DataProvider dataProvider;



        //<p>LabelProvider used for the Plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        LabelProvider labelProvider;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Plot(Plot & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~Plot(); 



        //<p>Method getting the unit from the data provider and setting the attribute unit with this value.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * getUnit(); 



        //<p>Method getting the title with the help of the data provider and the label provider and setting the title attribute.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * getTitle(); 



        //<p>Method that enables the rendering of the legend.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int showLegend(bool legendShow=false); 



        //<p>Method defining the legend with the help of the data provider and the label provider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        const char * getLegend(); 



        //<p>Method that will save the plotted plot in a file.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int write(const char * fileName); 



        //<p>Method that allows the user to set his own title.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int setTitle(const char * title); 



        //<p>Method allowing the user to set the PlottingStyle.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int setPlottingStyle(PlottingStyle style); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int setUnit(const char * unit); 



        //<p>Sets the data provider used for the plots.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int setDataProvider(DataProvider dataProvider); 



        //<p>Gets the data provider used.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        DataProvider getDataProvider(); 



        //<p>Sets the label provider used for the plots.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        int setLabelProvider(LabelProvider labelProvider); 



        //<p>Gets the label provider used.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        LabelProvider getLabelProvider(); 



};  //end class Plot



#endif
