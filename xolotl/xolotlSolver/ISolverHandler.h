#ifndef ISOLVERHANDLER_H
#define ISOLVERHANDLER_H

// Includes
#include <petscsys.h>
#include <petscts.h>
#include <petscdmda.h>
#include <memory>
#include <ITemperatureHandler.h>
#include <IDiffusionHandler.h>
#include <IAdvectionHandler.h>
#include <ITrapMutationHandler.h>
#include <IReSolutionHandler.h>
#include <IMaterialFactory.h>
#include <IReactionNetwork.h>

namespace xolotlSolver {

template<typename ValueType, typename SeedType>
class RandomNumberGenerator;

/**
 * Realizations of this interface are responsible for the actual implementation of
 * each piece of the solver. It is created to handle the multiple dimensions more easily.
 *
 * TODO - It will be better to have a PETSc-free solver handler interface (like ISolver)
 */
class ISolverHandler {

public:

	/**
	 * The destructor.
	 */
	virtual ~ISolverHandler() {
	}

	/**
	 * Initialize all the physics handlers that are needed to solve the ADR equations.
	 *
	 * @param material The material factory
	 * @param tempHandler The temperature handler
	 * @param options The Xolotl options
	 */
	virtual void initializeHandlers(
			std::shared_ptr<xolotlFactory::IMaterialFactory> material,
			std::shared_ptr<xolotlCore::ITemperatureHandler> tempHandler,
			const xolotlCore::Options &options) = 0;

	/**
	 * Create everything needed before starting to solve.
	 *
	 * @param da The PETSc distributed array
	 */
	virtual void createSolverContext(DM &da) = 0;

	/**
	 * Initialize the concentration solution vector.
	 *
	 * @param da The PETSc distributed array
	 * @param C The PETSc solution vector
	 */
	virtual void initializeConcentration(DM &da, Vec &C) = 0;

	/**
	 * Set the concentrations to 0.0 where the GBs are.
	 *
	 * @param da The PETSc distributed array
	 * @param C The PETSc solution vector
	 */
	virtual void initGBLocation(DM &da, Vec &C) = 0;

	/**
	 * Compute the new concentrations for the RHS function given an initial
	 * vector of concentrations.
	 *
	 * @param ts The PETSc time stepper
	 * @param localC The PETSc local solution vector
	 * @param F The updated PETSc solution vector
	 * @param ftime The real time
	 */
	virtual void updateConcentration(TS &ts, Vec &localC, Vec &F,
			PetscReal ftime) = 0;

	/**
	 * Compute the off-diagonal part of the Jacobian which is related to cluster's motion.
	 *
	 * @param ts The PETSc time stepper
	 * @param localC The PETSc local solution vector
	 * @param J The Jacobian
	 * @param ftime The real time
	 */
	virtual void computeOffDiagonalJacobian(TS &ts, Vec &localC, Mat &J,
			PetscReal ftime) = 0;

	/**
	 * Compute the diagonal part of the Jacobian which is related to cluster reactions.
	 *
	 * @param ts The PETSc time stepper
	 * @param localC The PETSc local solution vector
	 * @param J The Jacobian
	 * @param ftime The real time
	 */
	virtual void computeDiagonalJacobian(TS &ts, Vec &localC, Mat &J,
			PetscReal ftime) = 0;

	/**
	 * Get the grid in the x direction.
	 *
	 * @return The grid in the x direction
	 */
	virtual std::vector<double> getXGrid() const = 0;

	/**
	 * Get the step size in the y direction.
	 *
	 * @return The step size in the y direction
	 */
	virtual double getStepSizeY() const = 0;

	/**
	 * Get the step size in the z direction.
	 *
	 * @return The step size in the z direction
	 */
	virtual double getStepSizeZ() const = 0;

	/**
	 * Get the number of dimensions of the problem.
	 *
	 * @return The number of dimensions
	 */
	virtual int getDimension() const = 0;

	/**
	 * Get the position of the surface.
	 *
	 * @param j The index on the grid in the y direction
	 * @param k The index on the grid in the z direction
	 * @return The position of the surface at this y,z coordinates
	 */
	virtual int getSurfacePosition(int j = -1, int k = -1) const = 0;

	/**
	 * Set the position of the surface.
	 *
	 * @param pos The index of the position
	 * @param j The index on the grid in the y direction
	 * @param k The index on the grid in the z direction
	 */
	virtual void setSurfacePosition(int pos, int j = -1, int k = -1) = 0;

	/**
	 * Get the initial vacancy concentration.
	 *
	 * @return The initial vacancy concentration
	 */
	virtual double getInitialVConc() const = 0;

	/**
	 * Get the sputtering yield.
	 *
	 * @return The sputtering yield
	 */
	virtual double getSputteringYield() const = 0;

	/**
	 * Get the bursting depth parameter.
	 *
	 * @return The depth parameter
	 */
	virtual double getTauBursting() const = 0;

	/**
	 * Get the grid left offset.
	 *
	 * @return The offset
	 */
	virtual int getLeftOffset() const = 0;

	/**
	 * Get the grid right offset.
	 *
	 * @return The offset
	 */
	virtual int getRightOffset() const = 0;

	/**
	 * Create the local Xe rate vector.
	 *
	 * @param a The size in the x direction
	 * @param b The size in the y direction
	 * @param c The size in the y direction
	 */
	virtual void createLocalXeRate(int a, int b = 1, int c = 1) = 0;

	/**
	 * Set the latest value of the local Xe rate.
	 *
	 * @param rate The latest value of rate
	 * @param i The x coordinate of the location
	 * @param j The y coordinate of the location
	 * @param z The z coordinate of the location
	 */
	virtual void setLocalXeRate(double rate, int i, int j = 0, int k = 0) = 0;

	/**
	 * Get the value of the local Xe rate.
	 *
	 * @return The local vector of rates
	 */
	virtual std::vector<std::vector<std::vector<double> > > * getLocalXeRate() = 0;

	/**
	 * Create the local Xe conc vector.
	 *
	 * @param a The size in the x direction
	 * @param b The size in the y direction
	 * @param c The size in the y direction
	 */
	virtual void createLocalXeConc(int a, int b = 1, int c = 1) = 0;

	/**
	 * Set the latest value of the local Xe conc.
	 *
	 * @param rate The latest value of conc
	 * @param i The x coordinate of the location
	 * @param j The y coordinate of the location
	 * @param z The z coordinate of the location
	 */
	virtual void setLocalXeConc(double conc, int i, int j = 0, int k = 0) = 0;

	/**
	 * Get the value of the local Xe conc.
	 *
	 * @return The local vector of rates
	 */
	virtual std::vector<std::vector<std::vector<double> > > * getLocalXeConc() = 0;

	/**
	 * Set the coordinates covered by the local grid.
	 *
	 * @param xs, xm The start and width in the X direction on the local MPI process
	 * @param ys, ym The start and width in the Y direction on the local MPI process
	 * @param zs, zm The start and width in the Z direction on the local MPI process
	 */
	virtual void setLocalCoordinates(int xs, int xm, int ys = 0, int ym = 0, int zs = 0,
			int zm = 0) = 0;

	/**
	 * Get the coordinates covered by the local grid.
	 *
	 * @param xs, xm The start and width in the X direction on the local MPI process
	 * @param Mx The total width in the X direction
	 * @param ys, ym The start and width in the Y direction on the local MPI process
	 * @param My The total width in the Y direction
	 * @param zs, zm The start and width in the Z direction on the local MPI process
	 * @param Mz The total width in the Z direction
	 */
	virtual void getLocalCoordinates(int &xs, int &xm, int &Mx, int &ys,
			int &ym, int &My, int &zs, int &zm, int &Mz) = 0;

	/**
	 * To know if the surface should be able to move.
	 *
	 * @return True if the surface should be able to move.
	 */
	virtual bool moveSurface() const = 0;

	/**
	 * To know if the bubble bursting should be used.
	 *
	 * @return True if we want the bubble bursting.
	 */
	virtual bool burstBubbles() const = 0;

	/**
	 * Get the flux handler.
	 *
	 * @return The flux handler
	 */
	virtual xolotlCore::IFluxHandler *getFluxHandler() const = 0;

	/**
	 * Get the temperature handler.
	 *
	 * @return The temperature handler
	 */
	virtual xolotlCore::ITemperatureHandler *getTemperatureHandler() const = 0;

	/**
	 * Get the diffusion handler.
	 *
	 * @return The diffusion handler
	 */
	virtual xolotlCore::IDiffusionHandler *getDiffusionHandler() const = 0;

	/**
	 * Get the advection handler.
	 *
	 * @return The first advection handler
	 */
	virtual xolotlCore::IAdvectionHandler *getAdvectionHandler() const = 0;

	/**
	 * Get the advection handlers.
	 *
	 * @return The first advection handlers
	 */
	virtual std::vector<xolotlCore::IAdvectionHandler *> getAdvectionHandlers() const = 0;

	/**
	 * Get the modified trap-mutation handler.
	 *
	 * @return The modified trap-mutation handler
	 */
	virtual xolotlCore::ITrapMutationHandler *getMutationHandler() const = 0;

	/**
	 * Get the network.
	 *
	 * @return The network
	 */
	virtual xolotlCore::IReactionNetwork& getNetwork() const = 0;

	/**
	 * Get the network name.
	 *
	 * @return The network name
	 */
	virtual std::string getNetworkName() const = 0;

	/**
	 * Access the random number generator.
	 * The generator will have already been seeded.
	 *
	 * @return The RandomNumberGenerator object to use.
	 */
	virtual RandomNumberGenerator<int, unsigned int>& getRNG(void) const = 0;

	/**
	 * Set the file name containing the location of GB.
	 *
	 * @param name The filename
	 */
	virtual void setGBFileName(std::string name) = 0;

	/**
	 * Get the vector containing the location of GB.
	 *
	 * @return The GB vector
	 */
	virtual std::vector<std::tuple<int, int, int> > getGBVector() const = 0;

	/**
	 * Set the location of one GB grid point.
	 *
	 * @param i, j, k The coordinate of the GB
	 */
	virtual void setGBLocation(int i, int j = 0, int k = 0) = 0;

	/**
	 * Reset the GB vector.
	 */
	virtual void resetGBVector() = 0;

};
//end class ISolverHandler

} /* namespace xolotlSolver */
#endif
