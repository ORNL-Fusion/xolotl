#ifndef IDESORPTIONHANDLER_H
#define IDESORPTIONHANDLER_H

// Includes
#include <PSICluster.h>
#include <IReactionNetwork.h>
#include <memory>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for all the physical parts
 * for the desorption of hydrogen clusters close to the surface.
 * The solver call these methods to handle the desorption.
 */
class IDesorptionHandler {

public:

	/**
	 * The destructor
	 */
	virtual ~IDesorptionHandler() {
	}

	/**
	 * The initialize method generates an empty vector for desorption.
	 *
	 * @param grid The grid on the x axis
	 * @param ny The number of grid points in the Y direction
	 * @param hy The step size in the Y direction
	 * @param nz The number of grid points in the Z direction
	 * @param hz The step size in the Z direction
	 */
	virtual void initialize(std::vector<double> grid, int ny = 0, double hy =
			0.0, int nz = 0, double hz = 0.0) = 0;

	/**
	 * This method defines which cluster desorbs at each grid point.
	 *
	 * @param surfacePos The index of the position of the surface
	 * @param network The network
	 * @param grid The grid on the x axis
	 */
	virtual void initializeIndex1D(int surfacePos,
			const IReactionNetwork& network, std::vector<double> grid) = 0;

	/**
	 * This method defines which cluster desorbs at each grid point.
	 *
	 * @param surfacePos The index of the position of the surface
	 * @param network The network
	 * @param grid The grid on the x axis
	 * @param ny The number of grid points in the Y direction
	 * @param hy The step size in the Y direction
	 */
	virtual void initializeIndex2D(std::vector<int> surfacePos,
			const IReactionNetwork& network, std::vector<double> grid, int ny,
			double hy) = 0;

	/**
	 * This method defines which cluster desorbs at each grid point.
	 *
	 * @param surfacePos The index of the position of the surface
	 * @param network The network
	 * @param grid The grid on the x axis
	 * @param ny The number of grid points in the Y direction
	 * @param hy The step size in the Y direction
	 * @param nz The number of grid points in the Z direction
	 * @param hz The step size in the Z direction
	 */
	virtual void initializeIndex3D(std::vector<std::vector<int> > surfacePos,
			const IReactionNetwork& network, std::vector<double> grid, int ny,
			double hy, int nz, double hz) = 0;

	/**
	 * Compute the flux due to the desorption for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the desorption is computed
	 * @param updatedConcOffset The pointer to the array of the concentration
	 * at the grid point where the desorption is computed used to find the
	 * next solution
	 * @param xi The index of the position on the grid in the depth direction
	 * @param xs The beginning of the grid on this process
	 * @param yj The index of the position on the grid in the Y direction
	 * @param zk The index of the position on the grid in the Z direction
	 */
	virtual void computeDesorption(double *concOffset,
			double *updatedConcOffset, int xi, int xs, int yj = 0,
			int zk = 0) = 0;

	/**
	 * Compute the partials due to the desorption for all the
	 * clusters given the position index xi. Returns the number of
	 * clusters that are undergoing desorption at this grid point.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * @param val The pointer to the array that will contain the values of
	 * partials for the desorption
	 * @param indices The pointer to the array that will contain the indices
	 * of the clusters
	 * @param xi The index of the position on the grid in the depth direction
	 * @param xs The beginning of the grid on this process
	 * @param yj The index of the position on the grid in the Y direction
	 * @param zk The index of the position on the grid in the Z direction
	 *
	 * @return The number of clusters that go through desorption
	 * at this grid point
	 */
	virtual int computePartialsForDesorption(double *val, int *indices, int xi,
			int xs, int yj = 0, int zk = 0) = 0;

	/**
	 * Get the total number of clusters in the network that can undergo desorption.
	 *
	 * @return The number of clusters
	 */
	virtual int getNumberOfDesorbing() const = 0;

};
//end class IDesorptionHandler

} /* namespace xolotlCore */
#endif
