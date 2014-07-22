#ifndef FEFITFLUXHANDLER_H
#define FEFITFLUXHANDLER_H

#include "FluxHandler.h"

namespace xolotlSolver{

/**
 * This class realizes the IFluxHandler interface to calculate the incident and outgoing fluxes.
 */
class FeFitFluxHandler: public FluxHandler {

public:

	/**
	 * The constructor
	 */
	FeFitFluxHandler();

	/**
	 * The Destructor
	 */
	~FeFitFluxHandler() { }

	/**
	 * Function to calculate and store the incident flux values at each grid point
	 * @param numGridpoints The total number of grid points that will be used
	 * @param step The step size between grid points
	 */
	void initializeFluxHandler(int numGridpoints, double step);

}; //end class FeFitFluxHandler

}

#endif
