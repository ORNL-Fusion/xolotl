#ifndef W221FITFLUXHANDLER_H
#define W221FITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (221) oriented tungsten material.
 */
class W221FitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (221).
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1) return 0.0;

		// Compute the fit
		double value = 2.31116232 + 9.06474711 * x - 8.70839454 * pow(x, 2)
		+ 3.56696873 * pow(x, 3) - 0.842027679 * pow(x, 4) + 0.124415681 * pow(x, 5)
		- 0.0116876851 * pow(x, 6) + 6.78081089e-04 * pow(x, 7) - 2.21502063e-05 * pow(x, 8)
		+ 3.11533589e-07 * pow(x, 9);

		return value;
	}

public:

	/**
	 * The constructor
	 */
	W221FitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~W221FitFluxHandler() {}

};
//end class W221FitFluxHandler

}

#endif
