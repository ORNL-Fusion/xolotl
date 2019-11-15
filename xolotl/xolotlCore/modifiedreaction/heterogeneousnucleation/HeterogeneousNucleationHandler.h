#ifndef HETEROGENEOUSNUCLEATIONHANDLER_H
#define HETEROGENEOUSNUCLEATIONHANDLER_H

// Includes
#include <IHeterogeneousNucleationHandler.h>

namespace xolotlCore {

/**
 * This class realizes IReSolutionHandler interface, responsible for the
 * heterogeneous nucleation of xenon clusters.
 */
class HeterogeneousNucleationHandler: public IHeterogeneousNucleationHandler {

private:

	//! The nucleation rate
	double nucleationRate;

	//! The fission yield
	double fissionYield;

public:

	/**
	 * The constructor
	 */
	HeterogeneousNucleationHandler() :
			nucleationRate(0.0), fissionYield(0.0) {
	}

	/**
	 * The destructor
	 */
	~HeterogeneousNucleationHandler() {
	}

	/**
	 * The initialize method has to add connectivity between the Xe clusters,
	 * It must also initialize the rates of the reaction.
	 *
	 * \see IHeterogeneousNucleationHandler.h
	 */
	virtual void initialize(const IReactionNetwork& network);

	/**
	 * This method update the rate for the heterogeneous nucleation if the fission rate
	 * changed, it should be called when the flux changes for instance.
	 *
	 * \see IHeterogeneousNucleationHandler.h
	 */
	virtual void updateHeterogeneousNucleationRate(double rate);

	/**
	 * This method updates the fission yield.
	 *
	 * \see IHeterogeneousNucleationHandler.h
	 */
	void setFissionYield(double yield);

	/**
	 * Compute the flux due to the heterogeneous nucleation for all the clusters,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * Xe_1 + Xe_1 --> Xe_2
	 *
	 * F(Xe_1) = -2F(Xe_2) = -2*rate if C(Xe_1) > 2*rate
	 * else F(Xe_1) = -2F[Xe_2] = -C(Xe_1)
	 *
	 * \see IHeterogeneousNucleationHandler.h
	 */
	virtual void computeHeterogeneousNucleation(const IReactionNetwork& network,
			double *concOffset, double *updatedConcOffset, int xi, int xs,
			int yj = 0, int zk = 0);

	/**
	 * Compute the partials due to the heterogeneous nucleation for all the
	 * clusters given the position index xi.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * Xe_1 + Xe_1 --> Xe_2
	 *
	 * dF(Xe_1)/dC_(Xe_1) = -2dF(Xe_2)/dC_(Xe_1) = 0 if C(Xe_1) > 2*rate
	 * else dF(Xe_1)/dC_(Xe_1) = -2dF(Xe_2)/dC_(Xe_1) = -1
	 *
	 * \see IHeterogeneousNucleationHandler.h
	 */
	virtual bool computePartialsForHeterogeneousNucleation(
			const IReactionNetwork& network, double *val, int *indices, int xi,
			int xs, int yj = 0, int zk = 0);

};
//end class HeterogeneousNucleationHandler

} /* namespace xolotlCore */
#endif
