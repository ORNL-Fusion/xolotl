#ifndef W311FITFLUXHANDLER_H
#define W311FITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (311) oriented tungsten material.
 */
class W311FitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (311).
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1) return 0.0;

		// Compute the fit
		double value = 4.57176909 + 7.61971222 * x - 8.33433933 * pow(x, 2)
		+ 3.22362434 * pow(x, 3) - 0.649190013 * pow(x, 4) + 0.0740133639 * pow(x, 5)
		- 4.70056679e-03 * pow(x, 6) + 1.47690370e-04 * pow(x, 7) - 1.56743267e-06 * pow(x, 8);

		return value;
	}

public:

	/**
	 * The constructor
	 */
	W311FitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~W311FitFluxHandler() {}

};
//end class W311FitFluxHandler

}

#endif
