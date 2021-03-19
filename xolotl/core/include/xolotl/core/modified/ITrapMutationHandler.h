#ifndef ITRAPMUTATIONHANDLER_H
#define ITRAPMUTATIONHANDLER_H

// Includes
#include <memory>

#include <xolotl/config.h>
#include <xolotl/core/advection/IAdvectionHandler.h>
#include <xolotl/core/network/IReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace modified
{
/**
 * Realizations of this interface are responsible for all the physical parts
 * for the modified trap-mutation of small helium clusters close to the surface.
 * The solver call these methods to handle the modified trap-mutation.
 */
class ITrapMutationHandler
{
public:
	/**
	 * The destructor
	 */
	virtual ~ITrapMutationHandler()
	{
	}

	/**
	 * The initialize method has to add connectivity between the He clusters and
	 * HeV clusters of same number of He, and I. It must also initialize the
	 * rates of the reactions and call initializeIndex to define which
	 * trap-mutation is allowed at each grid point.
	 *
	 * @param network The network
	 * @param dfill Connectivity map
	 * @param nx The number of grid points in the X direction
	 * @param ny The number of grid points in the Y direction
	 * @param nz The number of grid points in the Z direction
	 */
	virtual void
	initialize(network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& dfill, int nx, int ny = 0,
		int nz = 0) = 0;

	/**
	 * This method defines which trap-mutation is allowed at each grid point.
	 *
	 * @param surfacePos The index of the position of the surface
	 * @param network The network
	 * @param advectionHandlers The vector of advection handlers
	 * @param grid The grid on the x axis
	 * @param nx The number of grid points in the X direction
	 * @param xs The beginning of the grid on this process
	 */
	virtual void
	initializeIndex1D(int surfacePos, network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs) = 0;

	/**
	 * This method defines which trap-mutation is allowed at each grid point.
	 *
	 * @param surfacePos The index of the position of the surface
	 * @param network The network
	 * @param advectionHandlers The vector of advection handlers
	 * @param grid The grid on the x axis
	 * @param nx The number of grid points in the X direction
	 * @param xs The beginning of the grid on this process
	 * @param ny The number of grid points in the Y direction
	 * @param hy The step size in the Y direction
	 * @param ys The beginning of the grid on this process
	 */
	virtual void
	initializeIndex2D(std::vector<int> surfacePos,
		network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny, double hy,
		int ys) = 0;

	/**
	 * This method defines which trap-mutation is allowed at each grid point.
	 *
	 * @param surfacePos The index of the position of the surface
	 * @param network The network
	 * @param advectionHandlers The vector of advection handlers
	 * @param grid The grid on the x axis
	 * @param nx The number of grid points in the X direction
	 * @param xs The beginning of the grid on this process
	 * @param ny The number of grid points in the Y direction
	 * @param hy The step size in the Y direction
	 * @param ys The beginning of the grid on this process
	 * @param nz The number of grid points in the Z direction
	 * @param hz The step size in the Z direction
	 * @param zs The beginning of the grid on this process
	 */
	virtual void
	initializeIndex3D(std::vector<std::vector<int>> surfacePos,
		network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny, double hy, int ys,
		int nz, double hz, int zs) = 0;

	/**
	 * This method updates the rate for the modified trap-mutation if the rates
	 * changed in the network, it should be called when temperature changes
	 * for instance.
	 *
	 * @param rate The largest rate in the network
	 */
	virtual void
	updateTrapMutationRate(const double rate) = 0;

	/**
	 * This method sets the boolean to remember if we want attenuation or not.
	 *
	 * @param isAttenuation True if we want attenuation
	 */
	virtual void
	setAttenuation(bool isAttenuation) = 0;

	/**
	 * This method updates the rate that makes the modified trap-mutation
	 * inefficient with time, depending on the total helium concentration.
	 *
	 * @param conc The concentration of helium
	 */
	virtual void
	updateDisappearingRate(double conc) = 0;

	/**
	 * Compute the flux due to the modified trap-mutation for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the solver.
	 *
	 * @param network The network
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the trap-mutation is computed
	 * @param updatedConcOffset The pointer to the array of the concentration
	 * at the grid point where the trap-mutation is computed
	 * @param xi The index of the position on the grid in the depth direction
	 * @param yj The index of the position on the grid in the Y direction
	 * @param zk The index of the position on the grid in the Z direction
	 */
	virtual void
	computeTrapMutation(network::IReactionNetwork& network, double* concOffset,
		double* updatedConcOffset, int xi, int yj = 0, int zk = 0) = 0;

	/**
	 * Compute the partials due to the modified trap-mutation for all the
	 * clusters given the position index xi. Returns the number of helium
	 * clusters that are undergoing trap-mutation at this grid point.
	 * This method is called by the RHSJacobian from the solver.
	 *
	 * @param network The network
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the trap-mutation is computed
	 * @param val The pointer to the array that will contain the values of
	 * partials for the trap-mutation
	 * @param indices The pointer to the array that will contain the indices
	 * of the clusters
	 * @param xi The index of the position on the grid in the depth direction
	 * @param yj The index of the position on the grid in the Y direction
	 * @param zk The index of the position on the grid in the Z direction
	 *
	 * @return The number of helium clusters that go through modified
	 * trap-mutation at this grid point
	 */
	virtual int
	computePartialsForTrapMutation(network::IReactionNetwork& network,
		double* concOffset, double* val, IdType* indices, int xi, int yj = 0,
		int zk = 0) = 0;

	/**
	 * Get the total number of clusters in the network that can undergo trap
	 * mutation.
	 *
	 * @return The number of clusters
	 */
	virtual int
	getNumberOfMutating() const = 0;

	/**
	 * Enable reduced Jacobian method.
	 *
	 * @param The boolean
	 */
	virtual void
	setEnableReducedJacobian(bool reduced) = 0;
};
// end class ITrapMutationHandler

} /* namespace modified */
} /* namespace core */
} /* namespace xolotl */
#endif
