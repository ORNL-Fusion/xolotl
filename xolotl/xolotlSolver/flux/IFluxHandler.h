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

	/**
	 * This operation increments the Helium fluence at the current time step.
	 * @param dt			The length of the time step
	 * @param step			The grid step size
	 * @return				The value of the Helium fluence at the current time step
	 */
	virtual double incrementHeFluence(double dt, double step) = 0;

	/**
	 * This operation returns the Helium fluence
	 * @return	The Helium fluence at current time step
	 */
	virtual double getHeFluence() const = 0;

	/**
	 * This operation sets the maximum value of the Helium fluence.
	 * @param	The maximim Helium fluence value
	 */
	virtual void setMaxHeFluence(double fluence) = 0;

	/**
	 * This function returns the maximum value of the Helium fluence.
	 */
	virtual double getMaxHeFluence() const = 0;

	/**
	 * This operation gets whether or not the maximum Helium fluence will be used
	 * @return	True if program will use the max He fluence, and false if it won't
	 */
	virtual bool getUsingMaxHeFluence() = 0;

	/**
	 * This operation sets the factor to change the Helium flux.
	 * @param flux	Helium flux value
	 */
	virtual void setHeFlux(double flux) = 0;

	/**
	 * This operation gets the factor that changes the Helium flux.
	 * @return	Helium flux value
	 */
	virtual double getHeFlux() const = 0;

}; //end class IFluxHandler

}

#endif
