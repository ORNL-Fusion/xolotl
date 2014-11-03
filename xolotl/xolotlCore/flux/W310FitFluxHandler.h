#ifndef W310FITFLUXHANDLER_H
#define W310FITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident and outgoing fluxes.
 */
class W310FitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (310).
	 *
	 * @param x The position where to evaluate the fit
	 * @return the evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1) return 0.0;

		// Compute the fit
		double value = 6.35816188 + 5.82240276 * x - 8.93251652 * pow(x, 2)
		+ 4.4166042 * pow(x, 3) - 1.20697765 * pow(x, 4) + 0.206429504 * pow(x, 5)
		- 0.0227671179 * pow(x, 6) + 1.57937040e-03 * pow(x, 7) - 6.27860509e-05 * pow(x, 8)
		+ 1.09090863e-06 * pow(x, 9);

		return value;
	}

public:

	/**
	 * The constructor
	 */
	W310FitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~W310FitFluxHandler() {}

};
//end class W310FitFluxHandler

}

#endif
