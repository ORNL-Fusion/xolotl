#ifndef FITFLUXHANDLER_H
#define FITFLUXHANDLER_H

#include "IFluxHandler.h"

namespace xolotlSolver{

/**
 * This class realizes the IFluxHandler interface to calculate the incident and outgoing fluxes.
 */
class FitFluxHandler: public IFluxHandler {

public:

	/**
	 * The constructor
	 */
	FitFluxHandler();

	/**
	 * The Destructor
	 */
	~FitFluxHandler();

	/**
	 * This operation returns the incident flux for a specific cluster composition,
	 * position, and time.
	 * @param compositionVec  The composition of the cluster
	 * @param position        The position of the cluster
	 * @param currentTime     The time
	 * @return incidentFlux   The incident flux at the given position and time of the cluster with
	 * the specified composition
	 */
	virtual double getIncidentFlux(std::vector<int> compositionVec,
			std::vector<double> position, double currentTime);

	/**
	 * Given a specific concentration, position, and time, this operation sets the outgoing
	 * flux to the specified amount.
	 * @param composition  The composition of the cluster
	 * @param position     The position of the cluster
	 * @param time         The time
	 * @return outgoingFlux  The outgoing flux at the given position and time of the cluster with
	 * the specified composition
	 */
	virtual void setOutgoingFlux(std::vector<int> compositionVec,
			std::vector<int> position, double time, double outgoingFlux);

}; //end class FitFluxHandler

}

#endif
