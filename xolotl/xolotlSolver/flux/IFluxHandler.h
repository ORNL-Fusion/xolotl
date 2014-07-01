#ifndef IFLUXHANDLER_H
#define IFLUXHANDLER_H

#include <vector>

namespace xolotlSolver{

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * and outgoing flux calculations.
 */
class IFluxHandler {

public:

	virtual ~IFluxHandler() { }

	/**
	 * Function to calculate and store the incident flux values at each grid point
	 * @param numGridpoints The total number of grid points that will be used
	 * @param step The step size between grid points
	 */
	virtual void initializeFluxHandler(int numGridpoints, double step) = 0;

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
			std::vector<double> position, double currentTime) = 0;

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
			std::vector<int> position, double time, double outgoingFlux) = 0;

}; //end class IFluxHandler

}

#endif
