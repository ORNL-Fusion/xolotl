#ifndef IFLUXHANDLER_H
#define IFLUXHANDLER_H

#include <vector>
#include <string>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * flux calculations.
 */
class IFluxHandler {

public:

	virtual ~IFluxHandler() { }

	/**
	 * Compute and store the incident flux values at each grid point.
	 *
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
	 *
	 * @param fileName The name of the file where the values are stored
	 */
	virtual void initializeTimeProfile(const std::string& fileName) = 0;

	/**
	 * This operation returns the incident flux vector.
	 *
	 * @param currentTime The time
	 * @return The incident flux vector
	 */
	virtual std::vector<double> getIncidentFluxVec(double currentTime) = 0;

	/**
	 * This operation increments the helium fluence at the current time step.
	 *
	 * @param dt The length of the time step
	 */
	virtual void incrementHeFluence(double dt) = 0;

	/**
	 * This operation returns the helium fluence.
	 *
	 * @return The helium fluence
	 */
	virtual double getHeFluence() const = 0;

	/**
	 * This operation sets the factor to change the intensity of the helium flux.
	 *
	 * @param flux The helium flux intensity
	 */
	virtual void setHeFlux(double flux) = 0;

	/**
	 * This operation gets the factor that changes the helium flux intensity.
	 *
	 * @return The helium flux value
	 */
	virtual double getHeFlux() const = 0;

}; //end class IFluxHandler

}

#endif
