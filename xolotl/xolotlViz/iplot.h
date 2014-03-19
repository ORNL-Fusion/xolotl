#ifndef IPLOT_H
#define IPLOT_H
//Begin section for file iplot.h
//TODO: Add definitions that you want preserved
//End section for file iplot.h
#include "plottingstyle.h"
#include "dataprovider.h"
#include "labelprovider.h"



class IDataProvider; //Dependency Generated Source:IPlot Target:IDataProvider


class ILabelProvider; //Dependency Generated Source:IPlot Target:ILabelProvider





//<p>IPlot describe the structure needed to be able to plot data provided by IDataProvider. The user interacts with it through different method where he/she could set the data to plot, title, legend, plotting style, etc.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class IPlot
{

    //Begin section for IPlot
    //TODO: Add attributes that you want preserved
    //End section for IPlot

    private:


        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        IDataProvider * iDataProvider;



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ILabelProvider * iLabelProvider;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual ~IPlot(); 



        //<p>Method that enables the rendering of the legend.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int showLegend(bool legendShow=false) = 0; 



        //<p>Method managing everything that is related to the rendering of a plot.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int render(bool isCumulative=false) = 0; 



        //<p>Method that will save the plotted plot in a file.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int write(const char * fileName) = 0; 



        //<p>Method that allows the user to set his own title.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setTitle(const char * title) = 0; 



        //<p>Method allowing the user to set the PlottingStyle.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setPlottingStyle(PlottingStyle style) = 0; 



        //<p>Method allowing the user to set the unit of the value plotted.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setUnit(const char * unit) = 0; 



        //<p>Sets the data provider used for the plots.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setDataProvider(DataProvider dataProvider) = 0; 



        //<p>Gets the data provider used.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual DataProvider getDataProvider() = 0; 



        //<p>Sets the label provider used for the plots.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual int setLabelProvider(LabelProvider labelProvider) = 0; 



        //<p>Gets the label provider used.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual LabelProvider getLabelProvider() = 0; 



        //<p>Method getting the unit from the data provider and setting the attribute unit with this value.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getUnit() = 0; 



        //<p>Method getting the title with the help of the data provider and the label provider and setting the title attribute.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getTitle() = 0; 



        //<p>Method defining the legend with the help of the data provider and the label provider.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        virtual const char * getLegend() = 0; 



};  //end class IPlot



#endif
