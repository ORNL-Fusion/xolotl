#ifndef DIFFUSION1DHANDLER_H
#define DIFFUSION1DHANDLER_H

// Includes
#include "DiffusionHandler.h"

namespace xolotlCore {

/**
 * This class is a subclass of DiffusionHandler for the isotropic diffusion of
 * mobile clusters in 1D.
 */
class Diffusion1DHandler: public DiffusionHandler {
private:

	//! The vector to know which clusters are diffusing where
	std::vector<std::vector<bool> > diffusionGrid;

public:

	//! The Constructor
	Diffusion1DHandler() {
	}

	//! The Destructor
	~Diffusion1DHandler() {
	}

	/**
	 * Initialize an vector of the dimension of the physical domain times the number of diffusion
	 * clusters. For each location, True means the cluster is diffusion, False means it is not.
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
	void initializeDiffusionGrid(
			std::vector<IAdvectionHandler *> advectionHandlers,
			std::vector<double> grid, int nx, int xs, int ny = 0, double hy =
					0.0, int ys = 0, int nz = 0, double hz = 0.0, int zs = 0)
					override;

	/**
	 * Compute the flux due to the diffusion for all the cluster that are diffusing,
	 * given the space parameters.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * If D is the diffusion coefficient, and C_l, C_r, C_m the left, right, and
	 * middle concentration of this cluster, the value of the flux is:
	 *
	 * D * (2.0 / [a * (a + b)]) * (C_l + [a/b] * C_r - [1.0 + (a/b)] * C_m)
	 *
	 * @param network The network
	 * @param concVector The pointer to the pointer of arrays of concentration at middle,
	 * left, and right grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the diffusion is computed used to find the next solution
	 * @param hxLeft The step size on the left side of the point in the x direction (a)
	 * @param hxRight The step size on the right side of the point in the x direction (b)
	 * @param ix The position on the x grid
	 * @param sy The space parameter, depending on the grid step size in the y direction
	 * @param iy The position on the y grid
	 * @param sz The space parameter, depending on the grid step size in the z direction
	 * @param iz The position on the z grid
	 */
	void computeDiffusion(experimental::IReactionNetwork& network,
			double **concVector, double *updatedConcOffset, double hxLeft,
			double hxRight, int ix, double sy = 0.0, int iy = 0,
			double sz = 0.0, int iz = 0) const override;

	/**
	 * Compute the partials due to the diffusion of all the diffusing clusters given
	 * the space parameters.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * Using the same notation as for computeDiffusion, the partial derivative on the
	 * left grid point should be:
	 *
	 * D * (2.0 / [a * (a + b)])
	 *
	 * on the right grid point:
	 *
	 * D * (2.0 / [b * (a + b)])
	 *
	 * and on the middle grid point:
	 *
	 * - D * (2.0 / [a*b])
	 *
	 * @param network The network
	 * @param val The pointer to the array that will contain the values of partials
	 * for the diffusion
	 * @param indices The pointer to the array that will contain the indices of the
	 * diffusing clusters in the network
	 * @param hxLeft The step size on the left side of the point in the x direction (a)
	 * @param hxRight The step size on the right side of the point in the x direction (b)
	 * @param ix The position on the x grid
	 * @param sy The space parameter, depending on the grid step size in the y direction
	 * @param iy The position on the y grid
	 * @param sz The space parameter, depending on the grid step size in the z direction
	 * @param iz The position on the z grid
	 */
	void computePartialsForDiffusion(
			experimental::IReactionNetwork& network, double *val,
			int *indices, double hxLeft, double hxRight, int ix,
			double sy = 0.0, int iy = 0, double sz = 0.0, int iz = 0) const
					override;

};
//end class Diffusion1DHandler

} /* end namespace xolotlCore */
#endif
