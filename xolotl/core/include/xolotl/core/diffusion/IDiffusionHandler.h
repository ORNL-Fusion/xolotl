#pragma once

// Includes
#include <memory>

#include <Kokkos_Core.hpp>

#include <xolotl/config.h>
#include <xolotl/core/Types.h>
#include <xolotl/core/advection/IAdvectionHandler.h>
#include <xolotl/core/network/IReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace diffusion
{
/**
 * Realizations of this interface are responsible for all the physical parts
 * for the isotropic diffusion of mobile cluster. The solver call these methods
 * to handle the diffusion.
 */
class IDiffusionHandler
{
public:
	/**
	 * The destructor
	 */
	virtual ~IDiffusionHandler()
	{
	}

	/**
	 * Initialize the off-diagonal part of the Jacobian. If this step is skipped
	 * it won't be possible to set the partial derivatives for the diffusion.
	 *
	 * @param network The network
	 * @param ofill Map of connectivity for diffusing clusters.
	 */
	virtual void
	initialize(network::IReactionNetwork& network,
		std::vector<RowColPair>& idPairs) = 0;

	/**
	 * Initialize an array of the dimension of the physical domain times the
	 * number of diffusing clusters. For each location, True means the cluster
	 * is actually diffusing, False means it is not. Clusters do not diffuse
	 * on advection sinks.
	 *
	 * @param advectionHandlers The vector of advection handlers
	 * @param grid The spatial grid in the depth direction
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
	initializeDiffusionGrid(
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny = 0, double hy = 0.0,
		int ys = 0, int nz = 0, double hz = 0.0, int zs = 0) = 0;

	/**
	 * Compute the flux due to the diffusion for all the clusters that are
	 * diffusing, given the space parameters. This method is called by the
	 * RHSFunction from the solver.
	 *
	 * @param network The network
	 * @param concVector The pointer to the pointer of arrays of concentration
	 * at middle/ left/right/bottom/top/front/back grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at
	 * the grid point where the diffusion is computed
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
	computeDiffusion(network::IReactionNetwork& network,
		const StencilConcArray& concVector,
		Kokkos::View<double*> updatedConcOffset, double hxLeft, double hxRight,
		int ix, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) const = 0;

	/**
	 * Compute the partials due to the diffusion of all the diffusing clusters
	 * given the space parameters. This method is called by the RHSJacobian from
	 * the solver.
	 *
	 * @param network The network
	 * @param val The pointer to the array that will contain the values of
	 * partials for the diffusion
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
	computePartialsForDiffusion(network::IReactionNetwork& network,
		Kokkos::View<double*> val, double hxLeft, double hxRight, int ix,
		double sy = 0.0, int iy = 0, double sz = 0.0, int iz = 0) const = 0;

	/**
	 * Get the total number of diffusing clusters in the network.
	 *
	 * @return The number of diffusing clusters
	 */
	virtual int
	getNumberOfDiffusing() const = 0;

	/**
	 * Get the vector of IDs of diffusing clusters in the network.
	 *
	 * @return The vector
	 */
	virtual std::vector<IdType>
	getDiffusingIds() const = 0;
};
// end class IDiffusionHandler

} /* end namespace diffusion */
} /* end namespace core */
} /* end namespace xolotl */
