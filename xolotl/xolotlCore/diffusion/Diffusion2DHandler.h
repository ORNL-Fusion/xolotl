#ifndef DIFFUSION2DHANDLER_H
#define DIFFUSION2DHANDLER_H

// Includes
#include "DiffusionHandler.h"

namespace xolotlCore {

/**
 * This class is a subclass of DiffusionHandler for the isotropic diffusion in 2D.
 */
class Diffusion2DHandler: public DiffusionHandler {
public:

	//! The Constructor
	Diffusion2DHandler() {}

	//! The Destructor
	~Diffusion2DHandler() {}

	/**
	 * Compute the flux due to the diffusion for all the cluster that are diffusing,
	 * given the space parameter s.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * If D is the diffusion coefficient, and C_l, C_r, C_b, C_t, C_m the left, right,
	 * bottom, top, and middle concentration of this cluster, the value to add to the
	 * updated concentration is:
	 *
	 * D * [sx * (C_l + C_r  - 4*C_m) + sy * (C_b + C_t - 2*C_m)]
	 *
	 * @param network The network
	 * @param concVector The pointer to the pointer of arrays of concentration at middle, left,
	 * right, bottom, and top grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the diffusion is computed used to find the next solution
	 * @param sx The space parameter, depending on the grid step size in the x direction
	 * @param sy The space parameter, depending on the grid step size in the y direction
	 * @param sz The space parameter, depending on the grid step size in the z direction
	 */
	void computeDiffusion(PSIClusterReactionNetwork *network,
			double **concVector, double *updatedConcOffset,
			double sx, double sy = 0.0, double sz = 0.0);

	/**
	 * Compute the partials due to the diffusion of all the diffusing clusters given
	 * the space parameter s.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * On the directly neighboring grid points the partial will be D * si with D the diffusion
	 * coefficient and the si the associated space parameter, and on this grid point the value
	 * of the partial will be -2 * D * (sx + sy).
	 *
	 * @param network The network
	 * @param val The pointer to the array that will contain the values of partials
	 * for the diffusion
	 * @param indices The pointer to the array that will contain the indices of the
	 * diffusing clusters in the network
	 * @param sx The space parameter, depending on the grid step size in the x direction
	 * @param sy The space parameter, depending on the grid step size in the y direction
	 * @param sz The space parameter, depending on the grid step size in the z direction
	 */
	void computePartialsForDiffusion(PSIClusterReactionNetwork *network,
			double *val, int *indices,
			double sx, double sy = 0.0, double sz = 0.0);

};
//end class Diffusion2DHandler

} /* end namespace xolotlCore */
#endif
