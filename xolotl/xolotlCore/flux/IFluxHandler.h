#ifndef IFLUXHANDLER_H
#define IFLUXHANDLER_H

#include <vector>
#include <string>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * and outgoing flux calculations.
 */
class IFluxHandler {

public:

	virtual ~IFluxHandler() { }

	/**
	 * Function to calculate and store the incident flux values at each grid point
	 * @param nx The total number of grid points that will be used on the x axis
	 * @param hx The step size between grid points on the x axis
	 * @param hy The step size between grid points on the y axis
	 * @param hz The step size between grid points on the z axis
	 */
	virtual void initializeFluxHandler(int nx, double hx, double hy = 1.0,
			double hz = 1.0) = 0;
	/**
	 * This method reads the values on the time profile file and store them in the
	 * time and amplitude vectors.
	 * @param fileName The name of the file where the values are stored
	 */
	virtual void initializeTimeProfile(std::string fileName) = 0;

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
	 * This operation returns the incident flux vector
	 * @param currentTime     	 The time
	 * @return incidentFluxVec   The incident flux vector
	 */
	virtual std::vector<double> getIncidentFluxVec(double currentTime) = 0;

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
	 */
	virtual void incrementHeFluence(double dt) = 0;

	/**
	 * This operation returns the Helium fluence
	 * @return	The Helium fluence at current time step
	 */
	virtual double getHeFluence() const = 0;

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
