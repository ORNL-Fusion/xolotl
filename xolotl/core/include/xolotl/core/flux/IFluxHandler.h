#pragma once

#include <string>
#include <vector>

#include <Kokkos_Core.hpp>

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
	 * the grid point where the diffusion is computed
	 * @param ix The position on the x grid
	 * @param surfacePos The current position of the surface
	 */
	virtual void
	computeIncidentFlux(double currentTime,
		Kokkos::View<double*> updatedConcOffset, int xi, int surfacePos) = 0;

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
	 * This operation sets the fluence.
	 *
	 * @param fluence The current fluence
	 */
	virtual void
	setFluence(std::vector<double> fluence) = 0;

	/**
	 * This operation returns the total fluence and effective fluences.
	 *
	 * @return The fluence
	 */
	virtual std::vector<double>
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

	/**
	 * Get the implanted flux for a specific sub network.
	 *
	 * @param The map of indices for this sub network
	 * @return The vector flux, first is the ID and second is the value.
	 */
	virtual std::vector<std::pair<IdType, double>>
	getImplantedFlux(std::vector<IdType> map) = 0;

	/**
	 * Set the implanted flux for each sub network.
	 *
	 * @param fluxVector With first is the ID and second is the value
	 */
	virtual void
	setImplantedFlux(std::vector<std::pair<IdType, double>> fluxVector) = 0;

	/**
	 * This operation gets the vector of flux amplitudes for each cluster at
	 * this time.
	 *
	 * @param time The current time
	 * @return The flux amplitudes
	 */
	virtual std::vector<double>
	getInstantFlux(double time) const = 0;

	/**
	 * This operation gets the vector of cluster IDs for each generated cluster.
	 *
	 * @return The flux indices
	 */
	virtual std::vector<IdType>
	getFluxIndices() const = 0;

	/**
	 * This operation gets the reduction factors each generated cluster.
	 *
	 * @return The factors
	 */
	virtual std::vector<double>
	getReductionFactors() const = 0;
};
// end class IFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl
