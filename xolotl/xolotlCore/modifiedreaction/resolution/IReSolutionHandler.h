#ifndef IRESOLUTIONHANDLER_H
#define IRESOLUTIONHANDLER_H

// Includes
#include <IReactionNetwork.h>
#include <memory>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for all the physical parts
 * for the re-solution of xenon clusterse.
 * The solver call these methods to handle the re-solution.
 */
class IReSolutionHandler {

public:

	/**
	 * The destructor
	 */
	virtual ~IReSolutionHandler() {
	}

	/**
	 * The initialize method has to add connectivity between the Xe clusters,
	 * It must also initialize the rates of the reactions.
	 *
	 * @param network The network
	 */
	virtual void initialize(const IReactionNetwork& network) = 0;

	/**
	 * This method update the rate for the re-solution if the fission rate
	 * changed, it should be called when the flux changes for instance.
	 *
	 * @param rate The fission rate
	 */
	virtual void updateReSolutionRate(double rate) = 0;

	/**
	 * Compute the flux due to the re-solution for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param network The network
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the trap-mutation is computed
	 * @param updatedConcOffset The pointer to the array of the concentration
	 * at the grid point where the trap-mutation is computed used to find the
	 * next solution
	 * @param xi The index of the position on the grid in the depth direction
	 * @param xs The beginning of the grid on this process
	 * @param yj The index of the position on the grid in the Y direction
	 * @param zk The index of the position on the grid in the Z direction
	 */
	virtual void computeReSolution(const IReactionNetwork& network,
			double *concOffset, double *updatedConcOffset, int xi, int xs,
			int yj = 0, int zk = 0) = 0;

	/**
	 * Compute the partials due to the re-solution for all the
	 * clusters given the position index xi. Returns the number of xenon
	 * clusters that are re-soluted at this grid point.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * @param network The network
	 * @param val The pointer to the array that will contain the values of
	 * partials for the trap-mutation
	 * @param indices The pointer to the array that will contain the indices
	 * of the clusters
	 * @param xi The index of the position on the grid in the depth direction
	 * @param xs The beginning of the grid on this process
	 * @param yj The index of the position on the grid in the Y direction
	 * @param zk The index of the position on the grid in the Z direction
	 *
	 * @return The number of xenon clusters that go re-soluted
	 * at this grid point
	 */
	virtual int computePartialsForReSolution(const IReactionNetwork& network,
			double *val, int *indices, int xi, int xs, int yj = 0,
			int zk = 0) = 0;

	/**
	 * Get the total number of clusters in the network that can undergo re-solution.
	 *
	 * @return The number of clusters
	 */
	virtual int getNumberOfReSoluting() const = 0;

};
//end class IReSolutionHandler

} /* namespace xolotlCore */
#endif
