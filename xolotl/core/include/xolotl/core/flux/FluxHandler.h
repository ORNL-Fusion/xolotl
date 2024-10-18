#pragma once

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
	 * View copy of incidentFluxVec
	 */
	Kokkos::View<double**> incidentFlux;

	/**
	 * The reduction factors for each deposition.
	 */
	std::vector<double> reductionFactors;

	/**
	 * Vector to hold the position at each grid
	 * point (x position).
	 */
	std::vector<double> xGrid;

	/**
	 *  Fluence.
	 */
	std::vector<double> fluence;

	/**
	 * The amplitude of the flux.
	 */
	double fluxAmplitude;

	/**
	 * The indices of the incoming clusters.
	 */
	std::vector<IdType> fluxIndices;

	/**
	 * View copy of fluxIndices
	 */
	Kokkos::View<IdType*> fluxIds;

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
	 * Value of the cascade dose.
	 */
	double cascadeDose;

	/**
	 * Value of remaining cascade efficiency.
	 */
	double cascadeEfficiency;

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

	/**
	 * This method copies flux indices to device view
	 */
	void
	syncFluxIndices();

	/**
	 * This method copies incident flux data to device view
	 */
	void
	syncIncidentFluxVec();

public:
	FluxHandler(const options::IOptions&);

	~FluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid) override;

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeTimeProfile(const std::string& fileName) final;

	/**
	 * \see IFluxHandler.h
	 */
	virtual void
	computeIncidentFlux(double currentTime, Kokkos::View<const double*>,
		Kokkos::View<double*> updatedConcOffset, int xi,
		int surfacePos) override;

	/**
	 * \see IFluxHandler.h
	 */
	void
	incrementFluence(double dt) override;

	/**
	 * \see IFluxHandler.h
	 */
	void
	computeFluence(double time) override;

	/**
	 * \see IFluxHandler.h
	 */
	void
	setFluence(std::vector<double> fluence) override;

	/**
	 * \see IFluxHandler.h
	 */
	std::vector<double>
	getFluence() const override;

	/**
	 * \see IFluxHandler.h
	 */
	void
	setFluxAmplitude(double flux) final;

	/**
	 * \see IFluxHandler.h
	 */
	double
	getFluxAmplitude() const override;

	/**
	 * \see IFluxHandler.h
	 */
	double
	getFluxRate() const override;

	/**
	 * \see IFluxHandler.h
	 */
	void
	setPulseTime(double time) override
	{
		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	setProportion(double a) override
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
	virtual std::vector<std::pair<IdType, double>>
	getImplantedFlux(std::vector<IdType> map)
	{
		return std::vector<std::pair<IdType, double>>();
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	setImplantedFlux(std::vector<std::pair<IdType, double>> fluxVector) override
	{
		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
	std::vector<double>
	getInstantFlux(double time) const override;

	/**
	 * \see IFluxHandler.h
	 */
	std::vector<IdType>
	getFluxIndices() const override;

	/**
	 * \see IFluxHandler.h
	 */
	std::vector<double>
	getReductionFactors() const override;
};
// end class FluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl
