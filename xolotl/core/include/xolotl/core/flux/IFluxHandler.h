#ifndef IFLUXHANDLER_H
#define IFLUXHANDLER_H

#include <string>
#include <vector>

#include <xolotl/core/network/IReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * Realizations of this interface are responsible for handling the incident
 * (incoming) flux calculations.
 */
class IFluxHandler
{
public:
	virtual ~IFluxHandler()
	{
	}

	/**
	 * Compute and store the incident flux values at each grid point.
	 *
	 * @param network The reaction network
	 * @param surfacePos The current position of the surface
	 * @param grid The grid on the x axis
	 */
	virtual void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid) = 0;

	/**
	 * This method reads the values on the time profile file and store them in
	 * the time and amplitude vectors.
	 *
	 * @param fileName The name of the file where the values are stored
	 */
	virtual void
	initializeTimeProfile(const std::string& fileName) = 0;

	/**
	 * This operation computes the flux due to incoming particles at a given
	 * grid point.
	 *
	 * @param currentTime The time
	 * @param updatedConcOffset The pointer to the array of the concentration at
	 * the grid point where the diffusion is computed used to find the next
	 * solution
	 * @param ix The position on the x grid
	 * @param surfacePos The current position of the surface
	 */
	virtual void
	computeIncidentFlux(double currentTime, double* updatedConcOffset, int xi,
		int surfacePos) = 0;

	/**
	 * This operation increments the fluence at the current time step.
	 *
	 * @param dt The length of the time step
	 */
	virtual void
	incrementFluence(double dt) = 0;

	/**
	 * This operation computes the fluence at the given time.
	 *
	 * @param time The current time
	 */
	virtual void
	computeFluence(double time) = 0;

	/**
	 * This operation returns the fluence.
	 *
	 * @return The fluence
	 */
	virtual double
	getFluence() const = 0;

	/**
	 * This operation sets the factor to change the intensity of the flux.
	 *
	 * @param flux The flux intensity
	 */
	virtual void
	setFluxAmplitude(double flux) = 0;

	/**
	 * This operation gets the factor that changes the flux
	 * intensity/amplitude.
	 *
	 * @return The flux amplitude
	 */
	virtual double
	getFluxAmplitude() const = 0;

	/**
	 * This operation gets the flux rate used for re-solution.
	 *
	 * @return The flux rate.
	 */
	virtual double
	getFluxRate() const = 0;

	/**
	 * This operation sets the time of the pulse.
	 *
	 * @param time The total time of one pulse
	 */
	virtual void
	setPulseTime(double time) = 0;

	/**
	 * This operation sets proportion of the pulse that is on.
	 *
	 * @param a The proprotion
	 */
	virtual void
	setProportion(double a) = 0;
};
// end class IFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
