#ifndef W210FITFLUXHANDLER_H
#define W210FITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlSolver {

/**
 * This class realizes the IFluxHandler interface to calculate the incident and outgoing fluxes.
 */
class W210FitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (210).
	 *
	 * @param x The position where to evaluate the fit
	 * @return the evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1) return 0.0;

		// Compute the fit
		double value = 6.12789477 + 6.95072358 * x - 10.1148409 * pow(x, 2)
		+ 4.90422349 * pow(x, 3) - 1.28776241 * pow(x, 4) + 0.205459896 * pow(x, 5)
		- 2.04925417e-02 * pow(x, 6) + 1.24916704e-03 * pow(x, 7) - 4.25758944e-05 * pow(x, 8)
		+ 6.21713210e-07 * pow(x, 9);

		return value;
	}

public:

	/**
	 * The constructor
	 */
	W210FitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~W210FitFluxHandler() {}

};
//end class W210FitFluxHandler

}

#endif
