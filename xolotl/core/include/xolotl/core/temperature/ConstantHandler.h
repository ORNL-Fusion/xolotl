#pragma once

#include <xolotl/core/temperature/TemperatureHandler.h>
#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of a temperature constant with time.
 */
class ConstantHandler : public TemperatureHandler
{
private:
	/**
	 * The temperature in Kelvin
	 */
	double temperature;

public:
	/**
	 * Construct with provided constant temperature
	 *
	 * @param constTemperature the temperature
	 */
	ConstantHandler(double constTemperature);

	/**
	 * Construct with options
	 */
	ConstantHandler(const options::IOptions& options);

	/**
	 * The destructor.
	 */
	virtual ~ConstantHandler();

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 * Here it is a constant temperature.
	 *
	 * \see ITemperatureHandler.h
	 */
	double
	getTemperature(const plsm::SpaceVector<double, 3>&, double) const override
	{
		return temperature;
	}

	/**
	 * This operation sets the temperature given by the solver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	void
	setTemperature(double* solution) override
	{
		return;
	}

	/**
	 * This operation sets the heat coefficient to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	void
	setHeatCoefficient(double coef) override
	{
		return;
	}

	/**
	 * This operation sets the heat conductivity to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	void
	setHeatConductivity(double cond) override
	{
		return;
	}

	/**
	 * This operation sets the surface position.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	void
	updateSurfacePosition(int surfacePos) override
	{
		return;
	}

	/**
	 * Compute the flux due to the heat equation.
	 * This method is called by the RHSFunction from the PetscSolver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	void
	computeTemperature(double** concVector, double* updatedConcOffset,
		double hxLeft, double hxRight, int xi, double sy = 0.0, int iy = 0,
		double sz = 0.0, int iz = 0) override
	{
		return;
	}

	/**
	 * Compute the partials due to the heat equation.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	bool
	computePartialsForTemperature(double* val, int* indices, double hxLeft,
		double hxRight, int xi, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) override
	{
		return false;
	}
};
} // namespace temperature
} // namespace core
} // namespace xolotl
