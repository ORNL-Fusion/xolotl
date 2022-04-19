#ifndef ITEMPERATUREHANDLER_H
#define ITEMPERATUREHANDLER_H

#include <memory>
#include <vector>

#include <Kokkos_View.hpp>

#include <plsm/SpaceVector.h>

#include <xolotl/config.h>
#include <xolotl/core/network/IReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
/**
 * Realizations of this interface are responsible for handling the temperature.
 */
class ITemperatureHandler
{
public:
	/**
	 * The destructor.
	 */
	virtual ~ITemperatureHandler()
	{
	}

	/**
	 * This operation initializes the variables that need to be
	 * depending on the type of handler used.
	 *
	 * @param dof The number of degrees of freedom
	 * @param ofillMap Map indicating row/column of diffusing variables in
	 * off-diagonal fill map.
	 * @param dfillMap Map indicating row/column of diffusing variables in
	 * diagonal fill map.
	 */
	virtual void
	initializeTemperature(const int dof,
		network::IReactionNetwork::SparseFillMap& ofillMap,
		network::IReactionNetwork::SparseFillMap& dfillMap) = 0;

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 *
	 * @param fraction The position fraction on the grid
	 * @param currentTime The time
	 * @return The temperature
	 */
	virtual double
	getTemperature(const plsm::SpaceVector<double, 3>& fraction,
		double currentTime) const = 0;

	/**
	 * This operation sets the temperature given by the solver.
	 *
	 * @param solution The pointer to the array of solutions
	 */
	virtual void
	setTemperature(double* solution) = 0;

	virtual void
	setTemperature(Kokkos::View<const double*> solution) = 0;

	/**
	 * This operation sets the heat coefficient to use in the equation.
	 *
	 * @param coef The heat coefficient
	 */
	virtual void
	setHeatCoefficient(double coef) = 0;

	/**
	 * This operation sets the heat conductivity to use in the equation.
	 *
	 * @param coef The heat conductivity
	 */
	virtual void
	setHeatConductivity(double cond) = 0;

	/**
	 * This operation sets the surface position.
	 *
	 * @param surfacePos The surface location
	 */
	virtual void
	updateSurfacePosition(int surfacePos) = 0;

	/**
	 * Compute the flux due to the heat equation.
	 * This method is called by the RHSFunction from the solver.
	 *
	 * @param concVector The pointer to the pointer of arrays of concentration
	 * at middle/ left/right grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at
	 * the grid point where the heat equation is computed
	 * @param hxLeft The step size on the left side of the point in the x
	 * direction
	 * @param hxRight The step size on the right side of the point in the x
	 * direction
	 * @param ix The position on the x grid
	 * @param sy The space parameter, depending on the grid step size in the y
	 * direction
	 * @param iy The position on the y grid
	 * @param sz The space parameter, depending on the grid step size in the z
	 * direction
	 * @param iz The position on the z grid
	 */
	virtual void
	computeTemperature(double** concVector, double* updatedConcOffset,
		double hxLeft, double hxRight, int xi, double sy = 0.0, int iy = 0,
		double sz = 0.0, int iz = 0) = 0;

	virtual void
	computeTemperature(Kokkos::View<const double*>* concVector,
		Kokkos::View<double*> updatedConcOffset, double hxLeft, double hxRight,
		int xi, double sy = 0.0, int iy = 0, double sz = 0.0, int iz = 0) = 0;

	/**
	 * Compute the partials due to the heat equation.
	 * This method is called by the RHSJacobian from the solver.
	 *
	 * @param val The pointer to the array that will contain the values of
	 * partials for the heat equation
	 * @param indices The pointer to the array that will contain the indices of
	 * the temperature in the network
	 * @param hxLeft The step size on the left side of the point in the x
	 * direction
	 * @param hxRight The step size on the right side of the point in the x
	 * direction
	 * @param ix The position on the x grid
	 * @param sy The space parameter, depending on the grid step size in the y
	 * direction
	 * @param iy The position on the y grid
	 * @param sz The space parameter, depending on the grid step size in the z
	 * direction
	 * @param iz The position on the z grid
	 * @return True if the partials were updated
	 */
	virtual bool
	computePartialsForTemperature(double* val, IdType* indices, double hxLeft,
		double hxRight, int xi, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) = 0;
};
// end class ITemperatureHandler

} // namespace temperature
} // namespace core
} // namespace xolotl

#endif
