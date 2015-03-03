#ifndef W320FITFLUXHANDLER_H
#define W320FITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (320) oriented tungsten material.
 */
class W320FitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (320).
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1) return 0.0;

		// Compute the fit
		double value = 5.12198605 + 6.61204266 * x - 7.4440834 * pow(x, 2)
		+ 2.79732401 * pow(x, 3) - 0.538911433 * pow(x, 4) + 0.0584367459 * pow(x, 5)
		- 3.53728409e-03 * pow(x, 6) + 1.07643241e-04 * pow(x, 7) - 1.17849962e-06 * pow(x, 8);

		return value;
	}

public:

	/**
	 * The constructor
	 */
	W320FitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~W320FitFluxHandler() {}

};
//end class W320FitFluxHandler

}

#endif
