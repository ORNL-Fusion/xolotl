#ifndef POINT_H
#define POINT_H
//Begin section for file point.h
//TODO: Add definitions that you want preserved
//End section for file point.h




//<p>Class describing the structure of data points. The attributes are the three spatial dimensions, the time, and the value of the quantity under consideration at this position.</p>
//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
class Point
{

    //Begin section for Point
    //TODO: Add attributes that you want preserved
    //End section for Point

    private:


        //<p>The time step.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double t;



        //<p>The X position on the grid.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double x;



        //<p>The Y position on the grid.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double y;



        //<p>The Z position on the grid.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double z;



        //<p>Value of the quantity of interest at the time step t and position on the grid x,y,z.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        double value;




    public:

        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        Point(Point & arg); 



        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        ~Point(); 



};  //end class Point



#endif
