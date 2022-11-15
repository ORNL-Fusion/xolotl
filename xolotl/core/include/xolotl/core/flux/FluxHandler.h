#ifndef FLUXHANDLER_H
#define FLUXHANDLER_H

#include <memory>
#include <vector>

#include <xolotl/core/Constants.h>
#include <xolotl/core/flux/IFluxHandler.h>
#include <xolotl/options/IOptions.h>

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
class FluxHandler : public IFluxHandler
{
protected:
	/**
	 * Vector to hold the incident flux values at each grid
	 * point (x position).
	 */
	std::vector<std::vector<double>> incidentFluxVec;

	/**
	 * Vector to hold the position at each grid
	 * point (x position).
	 */
	std::vector<double> xGrid;

	/**
	 *  Fluence.
	 */
	double fluence;

	/**
	 * The amplitude of the flux.
	 */
	double fluxAmplitude;

	/**
	 * The indices of the incoming clusters.
	 */
	std::vector<IdType> fluxIndices;

	/**
	 * Are we using a time profile for the amplitude of the incoming
	 * flux?
	 */
	bool useTimeProfile;

	/**
	 * Value of the fit function integrated on the grid.
	 */
	double normFactor;

	/**
	 * Vector to hold the time read from the input
	 * time profile file.
	 */
	std::vector<double> time;

	/**
	 * Vector to hold the amplitude read from the input
	 * time profile file.
	 */
	std::vector<double> amplitudes;

	/**
	 * Function that calculates the flux at a given position x (in nm).
	 * It needs to be implemented by the daughter classes.
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	virtual double
	FitFunction(double x)
	{
		return 0.0;
	}

	/**
	 * This method returns the value of the incident flux amplitude at
	 * the given time when a time profile is used.
	 *
	 * @param currentTime The time
	 * @return The value of the flux at this time
	 */
	double
	getProfileAmplitude(double currentTime) const;

	/**
	 * This method recomputes the values of the incident flux vector when
	 * conditions changed in the simulation.
	 *
	 * @param surfacePos The current position of the surface
	 */
	void
	recomputeFluxHandler(int surfacePos);

public:
	FluxHandler(const options::IOptions&);

	~FluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	virtual void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid);

	/**
	 * \see IFluxHandler.h
	 */
	virtual void
	initializeTimeProfile(const std::string& fileName);

	/**
	 * \see IFluxHandler.h
	 */
	virtual void
	computeIncidentFlux(double currentTime, double*, double* updatedConcOffset,
		int xi, int surfacePos);

	/**
	 * \see IFluxHandler.h
	 */
	virtual void
	incrementFluence(double dt);

	/**
	 * \see IFluxHandler.h
	 */
	virtual void
	computeFluence(double time);

	/**
	 * \see IFluxHandler.h
	 */
	virtual double
	getFluence() const;

	/**
	 * \see IFluxHandler.h
	 */
	virtual void
	setFluxAmplitude(double flux);

	/**
	 * \see IFluxHandler.h
	 */
	virtual double
	getFluxAmplitude() const;

	/**
	 * \see IFluxHandler.h
	 */
	virtual double
	getFluxRate() const;

	/**
	 * \see IFluxHandler.h
	 */
	virtual void
	setPulseTime(double time)
	{
		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
	virtual void
	setProportion(double a)
	{
		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
	virtual void
	setFissionYield(double yield)
	{
		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
	virtual std::vector<double>
	getInstantFlux(double time) const;

	/**
	 * \see IFluxHandler.h
	 */
	virtual std::vector<IdType>
	getFluxIndices() const;
};
// end class FluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
