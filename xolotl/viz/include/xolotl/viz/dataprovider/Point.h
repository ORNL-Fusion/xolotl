#ifndef POINT_H
#define POINT_H

// Includes
#include <string>

namespace xolotl {
namespace viz {
namespace dataprovider {

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
	double t = 0.;

	/**
	 * The X position on the grid.
	 */
	double x = 0.;

	/**
	 * The Y position on the grid.
	 */
	double y = 0.;

	/**
	 * The Z position on the grid.
	 */
	double z = 0.;

	/**
	 * Value of the quantity of interest at the time step t and position on the grid x,y,z.
	 */
	double value = 0.;

	/**
	 * The default constructor
	 */
	Point() {}

	/**
	 * The destructor
	 */
	~Point() {}

};

//end class Point

} /* namespace dataprovider */
} /* namespace viz */
} /* namespace xolotl */

#endif
