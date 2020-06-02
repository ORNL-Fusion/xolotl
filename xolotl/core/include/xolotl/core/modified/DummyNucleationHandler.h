#ifndef DUMMYNUCLEATIONHANDLER_H
#define DUMMYNUCLEATIONHANDLER_H

// Includes
#include <xolotl/core/modified/HeterogeneousNucleationHandler.h>

namespace xolotl {
namespace core {
namespace modified {

/**
 * This class realizes IReSolutionHandler interface, responsible for the
 * heterogeneous nucleation of xenon clusters. Here it is a dummy class,
 * meaning that it should not do anything.
 */
class DummyNucleationHandler: public HeterogeneousNucleationHandler {

public:

	/**
	 * The constructor
	 */
	DummyNucleationHandler() :
			HeterogeneousNucleationHandler() {
	}

	/**
	 * The destructor
	 */
	~DummyNucleationHandler() {
	}

	/**
	 * The initialize method has to add connectivity between the Xe clusters,
	 * It must also initialize the rates of the reaction.
	 *
	 * \see IHeterogeneousNucleationHandler.h
	 */
	void initialize(network::IReactionNetwork& network,
			network::IReactionNetwork::SparseFillMap& dfill) {
		// Doesn't do anything
		return;
	}

	/**
	 * This method update the rate for the heterogeneous nucleation if the fission rate
	 * changed, it should be called when the flux changes for instance.
	 *
	 * \see IHeterogeneousNucleationHandler.h
	 */
	void updateHeterogeneousNucleationRate(double rate) {
		// Doesn't do anything
		return;
	}

	/**
	 * This method updates the fission yield.
	 *
	 * \see IHeterogeneousNucleationHandler.h
	 */
	void setFissionYield(double yield) {
		return;
	}

	/**
	 * Compute the flux due to the heterogeneous nucleation for all the clusters,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * \see IHeterogeneousNucleationHandler.h
	 */
	void computeHeterogeneousNucleation(network::IReactionNetwork& network,
			double *concOffset, double *updatedConcOffset, int xi, int xs,
			int yj = 0, int zk = 0) {
		// Doesn't do anything
		return;
	}

	/**
	 * Compute the partials due to the heterogeneous nucleation for all the
	 * clusters given the position index xi.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * dF(Xe_1)/dC_(Xe_1) = -2dF(Xe_2)/dC_(Xe_1) = 0 if C(Xe_1) > 2*rate
	 * else dF(Xe_1)/dC_(Xe_1) = -2dF(Xe_2)/dC_(Xe_1) = -1
	 *
	 * \see IHeterogeneousNucleationHandler.h
	 */
	bool computePartialsForHeterogeneousNucleation(
			network::IReactionNetwork& network, double *concOffset,
			double *val, int *indices, int xi, int xs, int yj = 0, int zk = 0) {
		// Doesn't do anything
		return false;
	}

};
//end class DummyNucleationHandler

} /* namespace modified */
} /* namespace core */
} /* namespace xolotl */
#endif
