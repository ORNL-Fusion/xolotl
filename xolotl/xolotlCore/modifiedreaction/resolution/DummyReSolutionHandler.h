#ifndef DUMMYRESOLUTIONHANDLER_H
#define DUMMYRESOLUTIONHANDLER_H

// Includes
#include <ReSolutionHandler.h>
#include <Constants.h>

namespace xolotlCore {

/**
 * This class realizes IReSolutionHandler interface, responsible for the
 * re-solution of xenon clusters. Here it is a dummy class,
 * meaning that it should not do anything.
 */
class DummyReSolutionHandler: public ReSolutionHandler {

public:

	/**
	 * The constructor
	 */
	DummyReSolutionHandler() :
			ReSolutionHandler() {
	}

	/**
	 * The destructor
	 */
	~DummyReSolutionHandler() {
	}

	/**
	 * The initialize method has to add connectivity between the Xe clusters,
	 * It must also initialize the rates of the reactions.
	 *
	 * \see IReSolutionHandler.h
	 */
	void initialize(const IReactionNetwork& network,
			double electronicStoppingPower) {
		return;
	}

	/**
	 * This method updates the fission yield.
	 *
	 * \see IReSolutionHandler.h
	 */
	void setFissionYield(double yield) {
		return;
	}

	/**
	 * This method update the rate for the re-solution if the fission rate
	 * changed, it should be called when the flux changes for instance.
	 *
	 * \see IReSolutionHandler.h
	 */
	void updateReSolutionRate(double rate) {
		return;
	}

	/**
	 * Compute the flux due to the re-solution for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * Xe_i --> Xe_(i-n) + n Xe_1
	 *
	 * F(Xe_i) = -F[Xe_(i-n)] = -1/n * F(Xe_1) = -rate * C_(Xe_i)
	 *
	 * \see IReSolutionHandler.h
	 */
	void computeReSolution(const IReactionNetwork& network, double *concOffset,
			double *updatedConcOffset, int xi, int xs, int yj = 0, int zk = 0) {
		return;
	}

	/**
	 * Compute the partials due to the re-solution for all the
	 * clusters given the position index xi. Returns the number of xenon
	 * clusters that are re-soluted at this grid point.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * Xe_i --> Xe_(i-n) + n Xe_1
	 *
	 * dF(Xe_i)/dC_(Xe_i) = -dF[Xe_(i-n)]/dC_(Xe_i) = -1/n * dF(Xe_1)/dC_(Xe_i)
	 * 		= -rate
	 *
	 * \see IReSolutionHandler.h
	 */
	int computePartialsForReSolution(const IReactionNetwork& network,
			double *val, int *indices, int xi, int xs, int yj = 0, int zk = 0) {
		return 0;
	}

	/**
	 * Get the minimum size for a cluster to undergo re-solution.
	 *
	 * \see IReSolutionHandler.h
	 */
	int getMinSize() const {
		return 0;
	}

};
//end class DummyReSolutionHandler

} /* namespace xolotlCore */
#endif
