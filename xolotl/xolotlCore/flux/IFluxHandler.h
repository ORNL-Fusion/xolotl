#ifndef IFLUXHANDLER_H
#define IFLUXHANDLER_H

#include <vector>
#include <string>
#include <PSIClusterReactionNetwork.h>

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
	 * @param network The reaction network
	 * @param nx The total number of grid points that will be used on the x axis
	 * @param hx The step size between grid points on the x axis
	 */
	virtual void initializeFluxHandler(PSIClusterReactionNetwork *network,
			int nx, double hx) = 0;
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
	 * This operation returns the index of the cluster that is irradiating
	 * the material.
	 *
	 * @return The index of the incident flux cluster
	 */
	virtual int getIncidentFluxClusterIndex() = 0;

	/**
	 * This operation increments the fluence at the current time step.
	 *
	 * @param dt The length of the time step
	 */
	virtual void incrementFluence(double dt) = 0;

	/**
	 * This operation returns the fluence.
	 *
	 * @return The fluence
	 */
	virtual double getFluence() const = 0;

	/**
	 * This operation sets the factor to change the intensity of the flux.
	 *
	 * @param flux The flux intensity
	 */
	virtual void setFluxAmplitude(double flux) = 0;

	/**
	 * This operation gets the factor that changes the flux
	 * intensity/amplitude.
	 *
	 * @return The flux amplitude
	 */
	virtual double getFluxAmplitude() const = 0;

}; //end class IFluxHandler

}

#endif
