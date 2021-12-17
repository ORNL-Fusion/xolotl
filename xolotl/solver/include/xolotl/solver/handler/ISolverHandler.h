#ifndef ISOLVERHANDLER_H
#define ISOLVERHANDLER_H

// Includes
#include <petscdmda.h>
#include <petscsys.h>
#include <petscts.h>

#include <memory>

#include <xolotl/core/advection/IAdvectionHandler.h>
#include <xolotl/core/diffusion/IDiffusionHandler.h>
#include <xolotl/core/material/IMaterialHandler.h>
#include <xolotl/core/network/IReactionNetwork.h>
#include <xolotl/core/temperature/ITemperatureHandler.h>
#include <xolotl/options/IOptions.h>
#include <xolotl/perf/IPerfHandler.h>
#include <xolotl/util/Array.h>
#include <xolotl/util/RandomNumberGenerator.h>
#include <xolotl/viz/IVizHandler.h>

namespace xolotl
{
namespace solver
{
namespace handler
{
template <typename ValueType, typename SeedType>
class RandomNumberGenerator;

/**
 * Realizations of this interface are responsible for the actual implementation
 * of each piece of the solver. It is created to handle the multiple dimensions
 * more easily.
 *
 * TODO - It will be better to have a PETSc-free solver handler interface (like
 * ISolver)
 */
class ISolverHandler
{
public:
	/**
	 * The destructor.
	 */
	virtual ~ISolverHandler()
	{
	}

	/**
	 * Initialize all the physics handlers that are needed to solve the DR
	 * equations.
	 *
	 * @param material The material factory
	 * @param tempHandler The temperature handler
	 * @param options The Xolotl options
	 */
	virtual void
	initializeHandlers(core::material::IMaterialHandler* material,
		core::temperature::ITemperatureHandler* tempHandler,
		const options::IOptions& opts) = 0;

	/**
	 * Create everything needed before starting to solve.
	 *
	 * @param da The PETSc distributed array
	 */
	virtual void
	createSolverContext(DM& da) = 0;

	/**
	 * Initialize the concentration solution vector.
	 *
	 * @param da The PETSc distributed array
	 * @param C The PETSc solution vector
	 */
	virtual void
	initializeConcentration(DM& da, Vec& C) = 0;

	/**
	 * Set the concentrations to 0.0 where the GBs are.
	 *
	 * @param da The PETSc distributed array
	 * @param C The PETSc solution vector
	 */
	virtual void
	initGBLocation(DM& da, Vec& C) = 0;

	/**
	 * This operation get the concentration vector with the ids.
	 *
	 * @param da The PETSc distributed array
	 * @param C The PETSc solution vector
	 * @return The concentration vector
	 */
	virtual std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
	getConcVector(DM& da, Vec& C) = 0;

	/**
	 * This operation sets the concentration vector in the current state of the
	 * simulation.
	 *
	 * @param da The PETSc distributed array
	 * @param C The PETSc solution vector
	 * @param The concentration vector
	 */
	virtual void
	setConcVector(DM& da, Vec& C,
		std::vector<
			std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>&
			concVector) = 0;

	/**
	 * Get the previous time.
	 *
	 * @return The previous time
	 */
	virtual double
	getPreviousTime() = 0;

	/**
	 * Set the previous time.
	 *
	 * @param time The previous time
	 * @param updateFluence To know whether the fluence should be updated
	 */
	virtual void
	setPreviousTime(double time, bool updateFluence = false) = 0;

	/**
	 * Get the number of Xe that went to the GB.
	 *
	 * @return The number of Xenon
	 */
	virtual double
	getNXeGB() = 0;

	/**
	 * Set the number of Xe that went to the GB.
	 *
	 * @param nXe The number of Xenon
	 */
	virtual void
	setNXeGB(double nXe) = 0;

	/**
	 * Compute the new concentrations for the RHS function given an initial
	 * vector of concentrations.
	 *
	 * @param ts The PETSc time stepper
	 * @param localC The PETSc local solution vector
	 * @param F The updated PETSc solution vector
	 * @param ftime The real time
	 */
	virtual void
	updateConcentration(TS& ts, Vec& localC, Vec& F, PetscReal ftime) = 0;

	/**
	 * Compute the full Jacobian.
	 *
	 * @param ts The PETSc time stepper
	 * @param localC The PETSc local solution vector
	 * @param J The Jacobian
	 * @param ftime The real time
	 */
	virtual void
	computeJacobian(TS& ts, Vec& localC, Mat& J, PetscReal ftime) = 0;

	/**
	 * Get the grid in the x direction.
	 *
	 * @return The grid in the x direction
	 */
	virtual std::vector<double>
	getXGrid() const = 0;

	/**
	 * Get the step size in the y direction.
	 *
	 * @return The step size in the y direction
	 */
	virtual double
	getStepSizeY() const = 0;

	/**
	 * Get the step size in the z direction.
	 *
	 * @return The step size in the z direction
	 */
	virtual double
	getStepSizeZ() const = 0;

	/**
	 * Get the number of dimensions of the problem.
	 *
	 * @return The number of dimensions
	 */
	virtual int
	getDimension() const = 0;

	/**
	 * Get the position of the surface.
	 *
	 * @param j The index on the grid in the y direction
	 * @param k The index on the grid in the z direction
	 * @return The position of the surface at this y,z coordinates
	 */
	virtual IdType
	getSurfacePosition(IdType j = -1, IdType k = -1) const = 0;

	/**
	 * Set the position of the surface.
	 *
	 * @param pos The index of the position
	 * @param j The index on the grid in the y direction
	 * @param k The index on the grid in the z direction
	 */
	virtual void
	setSurfacePosition(IdType pos, IdType j = -1, IdType k = -1) = 0;

	/**
	 * Get the initial vacancy concentration.
	 *
	 * @return The initial vacancy concentration
	 */
	virtual double
	getInitialVConc() const = 0;

	/**
	 * Get the sputtering yield.
	 *
	 * @return The sputtering yield
	 */
	virtual double
	getSputteringYield() const = 0;

	/**
	 * Get the bursting depth parameter.
	 *
	 * @return The depth parameter
	 */
	virtual double
	getTauBursting() const = 0;

	/**
	 * Get the bursting factor for likelihood of bursting.
	 *
	 * @return The factor
	 */
	virtual double
	getBurstingFactor() const = 0;

	/**
	 * Get the HeV ratio.
	 *
	 * @return The ratio
	 */
	virtual double
	getHeVRatio() const = 0;

	/**
	 * Get the grid left offset.
	 *
	 * @return The offset
	 */
	virtual IdType
	getLeftOffset() const = 0;

	/**
	 * Get the grid right offset.
	 *
	 * @return The offset
	 */
	virtual IdType
	getRightOffset() const = 0;

	/**
	 * Create the local NE data vector.
	 *
	 * @param a The size in the x direction
	 * @param b The size in the y direction
	 * @param c The size in the y direction
	 */
	virtual void
	createLocalNE(IdType a, IdType b = 1, IdType c = 1) = 0;

	/**
	 * Set the latest value of the local Xe rate.
	 *
	 * @param rate The latest value of rate
	 * @param i The x coordinate of the location
	 * @param j The y coordinate of the location
	 * @param z The z coordinate of the location
	 */
	virtual void
	setLocalXeRate(double rate, IdType i, IdType j = 0, IdType k = 0) = 0;

	/**
	 * Set the whole vector of local NE data.
	 *
	 * @param rateVector The vector
	 */
	virtual void
	setLocalNE(
		const std::vector<std::vector<std::vector<std::array<double, 4>>>>&
			rateVector) = 0;

	/**
	 * Get the local NE data vector that needs to be passed to an app.
	 *
	 * @return The vector
	 */
	virtual std::vector<std::vector<std::vector<std::array<double, 4>>>>&
	getLocalNE() = 0;

	/**
	 * Set the latest value of the Xe flux.
	 *
	 * @param flux The latest value of flux
	 * @param i The x coordinate of the location
	 * @param j The y coordinate of the location
	 * @param z The z coordinate of the location
	 */
	virtual void
	setPreviousXeFlux(double flux, IdType i, IdType j = 0, IdType k = 0) = 0;

	/**
	 * Set the latest value of the Xe monomer concentration.
	 *
	 * @param conc The latest value of conc
	 * @param i The x coordinate of the location
	 * @param j The y coordinate of the location
	 * @param z The z coordinate of the location
	 */
	virtual void
	setMonomerConc(double conc, IdType i, IdType j = 0, IdType k = 0) = 0;

	/**
	 * Set the latest value of the volume fraction.
	 *
	 * @param frac The latest value of the fraction
	 * @param i The x coordinate of the location
	 * @param j The y coordinate of the location
	 * @param z The z coordinate of the location
	 */
	virtual void
	setVolumeFraction(double frac, IdType i, IdType j = 0, IdType k = 0) = 0;

	/**
	 * Set the coordinates covered by the local grid.
	 *
	 * @param xs, xm The start and width in the X direction on the local MPI
	 * process
	 * @param ys, ym The start and width in the Y direction on the local MPI
	 * process
	 * @param zs, zm The start and width in the Z direction on the local MPI
	 * process
	 */
	virtual void
	setLocalCoordinates(IdType xs, IdType xm, IdType ys = 0, IdType ym = 0,
		IdType zs = 0, IdType zm = 0) = 0;

	/**
	 * Get the coordinates covered by the local grid.
	 *
	 * @param xs, xm The start and width in the X direction on the local MPI
	 * process
	 * @param Mx The total width in the X direction
	 * @param ys, ym The start and width in the Y direction on the local MPI
	 * process
	 * @param My The total width in the Y direction
	 * @param zs, zm The start and width in the Z direction on the local MPI
	 * process
	 * @param Mz The total width in the Z direction
	 */
	virtual void
	getLocalCoordinates(IdType& xs, IdType& xm, IdType& Mx, IdType& ys,
		IdType& ym, IdType& My, IdType& zs, IdType& zm, IdType& Mz) = 0;

	/**
	 * To know if the surface should be able to move.
	 *
	 * @return True if the surface should be able to move.
	 */
	virtual bool
	moveSurface() const = 0;

	/**
	 * To know if the bubble bursting should be used.
	 *
	 * @return True if we want the bubble bursting.
	 */
	virtual bool
	burstBubbles() const = 0;

	/**
	 * To know if a temporal profile is used for the flux.
	 *
	 * @return True if temporal flux option is used.
	 */
	virtual bool
	temporalFlux() const = 0;

	/**
	 * Get the minimum size for computing average radius.
	 *
	 * @return The minimum size
	 */
	virtual std::vector<size_t>
	getMinSizes() const = 0;

	/**
	 * Get the flux handler.
	 *
	 * @return The flux handler
	 */
	virtual core::flux::IFluxHandler*
	getFluxHandler() const = 0;

	/**
	 * Get the temperature handler.
	 *
	 * @return The temperature handler
	 */
	virtual core::temperature::ITemperatureHandler*
	getTemperatureHandler() const = 0;

	/**
	 * Get the perf handler.
	 *
	 * @return The perf handler
	 */
	virtual std::shared_ptr<perf::IPerfHandler>
	getPerfHandler() const = 0;

	/**
	 * Get the viz handler.
	 *
	 * @return The viz handler
	 */
	virtual std::shared_ptr<viz::IVizHandler>
	getVizHandler() const = 0;

	/**
	 * Get the diffusion handler.
	 *
	 * @return The diffusion handler
	 */
	virtual core::diffusion::IDiffusionHandler*
	getDiffusionHandler() const = 0;

	/**
	 * Get the surface advection handler.
	 *
	 * @return The first advection handler
	 */
	virtual core::advection::IAdvectionHandler*
	getAdvectionHandler() const = 0;

	/**
	 * Get the advection handlers.
	 *
	 * @return The vector of handlers
	 */
	virtual std::vector<core::advection::IAdvectionHandler*>
	getAdvectionHandlers() const = 0;

	/**
	 * Get the network.
	 *
	 * @return The network
	 */
	virtual core::network::IReactionNetwork&
	getNetwork() const = 0;

	/**
	 * Get the network name.
	 *
	 * @return The network name
	 */
	virtual std::string
	getNetworkName() const = 0;

	/**
	 * Access the random number generator.
	 * The generator will have already been seeded.
	 *
	 * @return The RandomNumberGenerator object to use.
	 */
	virtual util::RandomNumberGenerator<int, unsigned int>&
	getRNG(void) const = 0;

	/**
	 * Get the vector containing the location of GB.
	 *
	 * @return The GB vector
	 */
	virtual std::vector<std::array<IdType, 3>>
	getGBVector() const = 0;

	/**
	 * Set the location of one GB grid point.
	 *
	 * @param i, j, k The coordinate of the GB
	 */
	virtual void
	setGBLocation(IdType i, IdType j = 0, IdType k = 0) = 0;

	/**
	 * Reset the GB vector.
	 */
	virtual void
	resetGBVector() = 0;
};
// end class ISolverHandler

} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */
#endif
