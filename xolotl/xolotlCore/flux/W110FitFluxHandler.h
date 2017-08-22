#ifndef W110FITFLUXHANDLER_H
#define W110FITFLUXHANDLER_H

#include "PSIFluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (110) oriented tungsten material.
 */
class W110FitFluxHandler: public PSIFluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (110).
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1) return 0.0;

		// Compute the fit
		double value = 7.93260868 + 1.49429886 * x - 4.48320209 * pow(x, 2)
		+ 1.97014869 * pow(x, 3) - 0.407986353 * pow(x, 4) + 0.0454535058 * pow(x, 5)
		- 0.0026618556 * pow(x, 6) + 0.0000678768532 * pow(x, 7) - 0.000000271171991 * pow(x, 8);

		return value;
	}

public:

	/**
	 * The constructor
	 */
	W110FitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~W110FitFluxHandler() {}

};
//end class W110FitFluxHandler

}

#endif
