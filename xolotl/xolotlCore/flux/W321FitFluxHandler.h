#ifndef W321FITFLUXHANDLER_H
#define W321FITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident and outgoing fluxes.
 */
class W321FitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (321).
	 *
	 * @param x The position where to evaluate the fit
	 * @return the evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1) return 0.0;

		// Compute the fit
		double value = 4.14077357 + 7.88582961 * x - 8.62196641 * pow(x, 2)
		+ 3.49053738 * pow(x, 3) - 0.758138247 * pow(x, 4) + 0.0966203298 * pow(x, 5)
		- 7.24058352e-03 * pow(x, 6) + 2.95631331e-04 * pow(x, 7) - 5.07710367e-06 * pow(x, 8);

		return value;
	}

public:

	/**
	 * The constructor
	 */
	W321FitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~W321FitFluxHandler() {}

};
//end class W321FitFluxHandler

}

#endif
