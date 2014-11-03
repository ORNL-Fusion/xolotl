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
	virtual void initializeOFill(std::shared_ptr<PSIClusterReactionNetwork> network,
			int *ofill) = 0;

	/**
	 * Compute the flux due to the diffusion for all the cluster that are diffusing,
	 * given the space parameter sx.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param network The network
	 * @param sx The space parameter, depending on the grid step size
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the diffusion is computed
	 * @param leftConcOffset The pointer to the array of concentration at the grid
	 * point to the left of where the diffusion is computed
	 * @param rightConcOffset The pointer to the array of concentration at the grid
	 * point to the right of where the diffusion is computed
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the diffusion is computed used to find the next solution
	 */
	virtual void computeDiffusion(std::shared_ptr<PSIClusterReactionNetwork> network, double sx,
			double *concOffset, double *leftConcOffset,
			double *rightConcOffset, double *updatedConcOffset) = 0;

	/**
	 * Compute the partials due to the diffusion of all the diffusing clusters given
	 * the space parameter sx.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * @param network The network
	 * @param sx The space parameter, depending on the grid step size
	 * @param val The pointer to the array that will contain the values of partials
	 * for the diffusion
	 * @param row The pointer to the array that will contain the indices of the row
	 * for the Jacobian
	 * @param col The pointer to the array that will contain the indices of the columns
	 * for the Jacobian
	 * @param xi The index of the grip point
	 * @param xs The index of the first grid point on the locally owned grid
	 */
	virtual void computePartialsForDiffusion(std::shared_ptr<PSIClusterReactionNetwork> network,
			double sx, double *val, int *row, int *col, int xi,
			int xs) = 0;

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
