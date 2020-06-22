#ifndef IADVECTIONHANDLER_H
#define IADVECTIONHANDLER_H

// Includes
#include <array>
#include <memory>

#include <plsm/SpaceVector.h>

#include <xolotl/core/network/IReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace advection
{
/**
 * Realizations of this interface are responsible for all the physical parts
 * for the advection of mobile helium cluster. The solver call these methods to
 * handle the advection.
 */
class IAdvectionHandler
{
public:
	/**
	 * The destructor
	 */
	virtual ~IAdvectionHandler()
	{
	}

	/**
	 * This function initialize the list of clusters that will move through
	 * advection and their corresponding sink strength (or driving forces).
	 *
	 * @param network The network
	 * @param ofillMap Map of connectivity for advecting clusters.
	 */
	virtual void
	initialize(network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& ofillMap) = 0;

	/**
	 * Set the number of dimension
	 *
	 * @param dim The number of dimensions of the problem
	 */
	virtual void
	setDimension(int dim) = 0;

	/**
	 * Set the location of the sink.
	 *
	 * @param pos The location of the sink
	 */
	virtual void
	setLocation(double pos) = 0;

	/**
	 * Initialize an array of the dimension of the physical domain times the
	 * number of advecting clusters. For each location, True means the cluster
	 * is moving, False means it is not.
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
	initializeAdvectionGrid(std::vector<IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny = 1, double hy = 0.0,
		int ys = 0, int nz = 1, double hz = 0.0, int zs = 0) = 0;

	/**
	 * Compute the flux due to the advection for all the helium clusters,
	 * given the space parameters and the position.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param network The network
	 * @param pos The position on the grid
	 * @param concVector The pointer to the pointer of arrays of concentration
	 * at middle, left, and right grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at
	 * the grid point where the advection is computed used to find the next
	 * solution
	 * @param hxLeft The step size on the left side of the point in the x
	 * direction
	 * @param hxRight The step size on the right side of the point in the x
	 * direction
	 * @param ix The position on the x grid
	 * @param hy The step size in the y direction
	 * @param iy The position on the y grid
	 * @param hz The step size in the z direction
	 * @param iz The position on the z grid
	 */
	virtual void
	computeAdvection(network::IReactionNetwork& network,
		const plsm::SpaceVector<double, 3>& pos, double** concVector,
		double* updatedConcOffset, double hxLeft, double hxRight, int ix,
		double hy = 0.0, int iy = 0, double hz = 0.0, int iz = 0) const = 0;

	/**
	 * Compute the partial derivatives due to the advection of all the helium
	 * clusters given the space parameters and the position. This method is
	 * called by the RHSJacobian from the PetscSolver.
	 *
	 * @param network The network
	 * @param h The space parameters in the three directions
	 * @param val The pointer to the array that will contain the values of
	 * partials for the advection
	 * @param indices The pointer to the array that will contain the indices of
	 * the advecting cluster in the network
	 * @param pos The position on the grid
	 * @param hxLeft The step size on the left side of the point in the x
	 * direction
	 * @param hxRight The step size on the right side of the point in the x
	 * direction
	 * @param ix The position on the x grid
	 * @param hy The step size in the y direction
	 * @param iy The position on the y grid
	 * @param hz The step size in the z direction
	 * @param iz The position on the z grid
	 */
	virtual void
	computePartialsForAdvection(network::IReactionNetwork& network, double* val,
		int* indices, const plsm::SpaceVector<double, 3>& pos, double hxLeft,
		double hxRight, int ix, double hy = 0.0, int iy = 0, double hz = 0.0,
		int iz = 0) const = 0;

	/**
	 * Compute the indices that will determine where the partial derivatives
	 * will be put in the Jacobian. This method is called by the RHSJacobian
	 * from the PetscSolver.
	 *
	 * @param pos The position on the grid
	 * @return The indices for the position in the Jacobian
	 */
	virtual std::array<int, 3>
	getStencilForAdvection(const plsm::SpaceVector<double, 3>& pos) const = 0;

	/**
	 * Check whether the grid point is located on the sink surface or not.
	 *
	 * @param pos The position on the grid
	 * @return True if the point is on the sink
	 */
	virtual bool
	isPointOnSink(const plsm::SpaceVector<double, 3>& pos) const = 0;

	/**
	 * Get the total number of advecting clusters in the network.
	 *
	 * @return The number of advecting clusters
	 */
	virtual int
	getNumberOfAdvecting() const = 0;

	/**
	 * Get the vector of index of advecting clusters in the network.
	 *
	 * @return The vector of advecting clusters
	 */
	virtual const std::vector<std::size_t>&
	getAdvectingClusters() = 0;

	/**
	 * Get the location of the sink.
	 *
	 * @return The location
	 */
	virtual double
	getLocation() const = 0;
};
// end class IAdvectionHandler

} /* namespace advection */
} /* namespace core */
} /* namespace xolotl */
#endif
