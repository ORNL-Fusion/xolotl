#ifndef IDIFFUSIONHANDLER_H
#define IDIFFUSIONHANDLER_H

// Includes
#include <PSICluster.h>
#include <PSIClusterReactionNetwork.h>
#include <memory>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for all the physical parts
 * for the diffusion of mobile cluster. The solver call these methods to handle
 * the diffusion.
 */
class IDiffusionHandler {

public:

	/**
	 * The destructor
	 */
	virtual ~IDiffusionHandler() {}

	/**
	 * Initialize the off-diagonal part of the Jacobian. If this step is skipped it
	 * won't be possible to set the partials for the diffusion.
	 *
	 * @param network The network
	 * @param ofill The pointer to the array that will contain the value 1 at the indices
	 * of the diffusing clusters, 0 if they are not diffusing
	 */
	virtual void initializeOFill(PSIClusterReactionNetwork *network,
			int *ofill) = 0;

	/**
	 * Compute the flux due to the diffusion for all the cluster that are diffusing,
	 * given the space parameter s.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param network The network
	 * @param concVector The pointer to the pointer of arrays of concentration at middle/
	 * left/right/bottom/top/front/back grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the diffusion is computed used to find the next solution
	 * @param sx The space parameter, depending on the grid step size in the x direction
	 * @param sy The space parameter, depending on the grid step size in the y direction
	 * @param sz The space parameter, depending on the grid step size in the z direction
	 */
	virtual void computeDiffusion(PSIClusterReactionNetwork *network,
			double **concVector, double *updatedConcOffset,
			double sx, double sy = 0.0, double sz = 0.0) = 0;

	/**
	 * Compute the partials due to the diffusion of all the diffusing clusters given
	 * the space parameter s.
	 * This method is called by the RHSJacobian from the PetscSolver.
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
	virtual void computePartialsForDiffusion(PSIClusterReactionNetwork *network,
			double *val, int *indices, double sx, double sy = 0.0, double sz = 0.0) = 0;

	/**
	 * Get the total number of diffusing clusters in the network.
	 *
	 * @return The number of diffusing clusters
	 */
	virtual int getNumberOfDiffusing() = 0;

};
//end class IDiffusionHandler

} /* namespace xolotlCore */
#endif
