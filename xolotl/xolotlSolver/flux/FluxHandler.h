#ifndef FLUXHANDLER_H
#define FLUXHANDLER_H

#include "IFluxHandler.h"
#include <vector>

namespace xolotlSolver{

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * and outgoing flux calculations.
 */
class FluxHandler: public IFluxHandler {

protected:

	/**
	 * Vector to hold the incident flux values at each grid
	 * point (x position)
	 */
	std::vector<double> incidentFluxVec;

	/**
	 * Step size between each grid point
	 */
	double stepSize;

public:

	FluxHandler()
		: stepSize(0.0e-16)
	{ }

	~FluxHandler() { }

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

}; //end class FluxHandler

}

#endif
