#ifndef TEMPERATUREGRADIENTHANDLER_H
#define TEMPERATUREGRADIENTHANDLER_H

#include <xolotl/core/temperature/ITemperatureHandler.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of a temperature constant with time but changing with location.
 */
class TemperatureGradientHandler : public ITemperatureHandler
{
private:
	/**
	 * The surface temperature in Kelvin
	 */
	double surfaceTemperature;

	/**
	 * The bulk temperature in Kelvin
	 */
	double bulkTemperature;

	/**
	 * The number of degrees of freedom in the network
	 */
	int dof;

	/**
	 * The default constructor is private because the TemperatureHandler
	 * must be initialized with a temperature
	 */
	TemperatureGradientHandler() :
		surfaceTemperature(0.0),
		bulkTemperature(0.0),
		dof(0)
	{
	}

public:
	/**
	 * The constructor
	 *
	 * @param temp The surface temperature
	 * @param grad The bulk temperature
	 */
	TemperatureGradientHandler(double surfTemp, double bulkTemp) :
		surfaceTemperature(surfTemp),
		bulkTemperature(bulkTemp),
		dof(0)
	{
	}

	/**
	 * The destructor.
	 */
	virtual ~TemperatureGradientHandler()
	{
	}

	/**
	 * This operation initializes the ofill and dfill arrays so that the
	 * temperature is connected correctly in the solver.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void
	initializeTemperature(const int _dof,
		network::IReactionNetwork::SparseFillMap& ofillMap,
		network::IReactionNetwork::SparseFillMap& dfillMap)
	{
		// Set dof
		dof = _dof;

		// Add the temperature to ofill
		ofillMap[dof].emplace_back(dof);

		// Add the temperature to dfill
		dfillMap[dof].emplace_back(dof);

		return;
	}

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 *
	 * @return The temperature
	 */
	virtual double
	getTemperature(const plsm::SpaceVector<double, 3>& fraction, double) const
	{
		return surfaceTemperature +
			(bulkTemperature - surfaceTemperature) * fraction[0];
	}

	/**
	 * This operation sets the temperature given by the solver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void
	setTemperature(double* solution)
	{
		return;
	}

	/**
	 * This operation sets the heat coefficient to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void
	setHeatCoefficient(double coef)
	{
		return;
	}

	/**
	 * This operation sets the heat conductivity to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void
	setHeatConductivity(double cond)
	{
		return;
	}

	/**
	 * This operation sets the surface position.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void
	updateSurfacePosition(int surfacePos)
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
	virtual void
	computeTemperature(double** concVector, double* updatedConcOffset,
		double hxLeft, double hxRight, int xi, double sy = 0.0, int iy = 0,
		double sz = 0.0, int iz = 0)
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
	virtual bool
	computePartialsForTemperature(double* val, int* indices, double hxLeft,
		double hxRight, int xi, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0)
	{
		return false;
	}
};
// end class TemperatureGradientHandler

} // namespace temperature
} // namespace core
} // namespace xolotl

#endif
