#ifndef POINT_H
#define POINT_H

// Includes
#include <string>

namespace xolotlViz {
namespace xolotl {
namespace viz {

/**
 * Class describing the structure of data points.
 * The attributes are the three spatial dimensions, the time, and the value of the quantity
 * under consideration at this position.
 */
class Point {

public:

	/**
	 * The time step.
	 */
	double t;

	/**
	 * The X position on the grid.
	 */
	double x;

	/**
	 * The Y position on the grid.
	 */
	double y;

	/**
	 * The Z position on the grid.
	 */
	double z;

	/**
	 * Value of the quantity of interest at the time step t and position on the grid x,y,z.
	 */
	double value;

	/**
	 * The default constructor
	 */
	Point() :
			t(0.0), x(0.0), y(0.0), z(0.0), value(0.0) {}

	/**
	 * The destructor
	 */
	~Point() {}

};

//end class Point

} /* namespace viz */
} /* namespace xolotl */

#endif
