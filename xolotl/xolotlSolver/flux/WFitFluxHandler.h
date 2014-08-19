#ifndef WFITFLUXHANDLER_H
#define WFITFLUXHANDLER_H

#include "FluxHandler.h"

namespace xolotlSolver {

/**
 * This class realizes the IFluxHandler interface to calculate the incident and outgoing fluxes.
 */
class WFitFluxHandler: public FluxHandler {
private:
	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (100).
	 *
	 * @param x The position where to evaluate he fit
	 * @return the evaluated value
	 */
	double FitFunction100(double x);

public:

	/**
	 * The constructor
	 */
	WFitFluxHandler() { }

	/**
	 * The Destructor
	 */
	~WFitFluxHandler() {
	}

	/**
	 * Function to calculate and store the incident flux values at each grid point
	 *
	 * @param numGridpoints The total number of grid points that will be used
	 * @param step The step size between grid points
	 */
	void initializeFluxHandler(int numGridpoints, double step);

};
//end class WFitFluxHandler

}

#endif
