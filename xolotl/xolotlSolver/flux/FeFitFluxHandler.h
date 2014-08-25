#ifndef FEFITFLUXHANDLER_H
#define FEFITFLUXHANDLER_H

#include "FluxHandler.h"

namespace xolotlSolver{

/**
 * This class realizes the IFluxHandler interface to calculate the incident and outgoing fluxes.
 */
class FeFitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized.
	 *
	 * @param x The position where to evaluate he fit
	 * @return the evaluated value
	 */
	double FitFunction(double x){
		// Value at which the flux goes to 0
		double x1 = 128.0;

		if (x > x1) return 0.0;

		// Set fit parameters
		double a0 = -0.00073090;
		double a1 = -0.0029330;
		double b1 = 0.0039810;
		double a2 = -0.0041960;
		double b2 = 0.0084180;
		double a3 = -0.0015640;
		double b3 = 0.0099430;
		double a4 = 0.0025910;
		double b4 = 0.0063010;
		double a5 = 0.0038630;
		double b5 = 0.000880;
		double a6 = 0.0022260;
		double b6 = -0.0017580;
		double a7 = 0.00053690;
		double b7 = -0.0013570;
		double a8 = 1.09200000;
		double b8 = -0.00035910;
		double w = 0.013880;

		// Compute the fit
		double value = a0 + a1 * cos(x * w) + b1 * sin(x * w)
		+ a2 * cos(2.0 * x * w) + b2 * sin(2.0 * x * w)
		+ a3 * cos(3.0 * x * w) + b3 * sin(3.0 * x * w)
		+ a4 * cos(4.0 * x * w) + b4 * sin(4.0 * x * w)
		+ a5 * cos(5.0 * x * w) + b5 * sin(5.0 * x * w)
		+ a6 * cos(6.0 * x * w) + b6 * sin(6.0 * x * w)
		+ a7 * cos(7.0 * x * w) + b7 * sin(7.0 * x * w)
		+ a8 * cos(8.0 * x * w) + b8 * sin(8.0 * x * w);

		return std::max(value, 0.0);
	}

public:

	/**
	 * The constructor
	 */
	FeFitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~FeFitFluxHandler() {}

}; //end class FeFitFluxHandler

}

#endif
