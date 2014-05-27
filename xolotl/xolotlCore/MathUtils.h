/*
 * MathUtils.h
 *
 * Various math utilities.
 *
 *  Created on: May 20, 2014
 *      Author: Jay Jay Billings
 */

#ifndef MATHUTILS_H_
#define MATHUTILS_H_

#include <limits>
#include <cmath>

namespace xolotlCore {

/**
 * This function tests two doubles to see if they are equal.
 * @param a The first double
 * @param b The secound double
 * @return True if the doubles are equal to within machine precision, false
 * otherwise.
 */
inline bool equal(double a, double b) {
	return std::fabs(b-a) < std::numeric_limits<double>::epsilon();
}

}

#endif /* MATHUTILS_H_ */
