#ifndef RESOLUTIONHANDLER_H
#define RESOLUTIONHANDLER_H

// Includes
#include <IReSolutionHandler.h>
#include <Constants.h>

namespace xolotlCore {

/**
 * This class realizes IReSolutionHandler interface, responsible for the
 * re-solution of xenon clusters.
 */
class ReSolutionHandler: public IReSolutionHandler {
protected:

	struct ReSolutionBase {

		/**
		 * The larger cluster in the pair
		 */
		IReactant* larger;

		/**
		 * The smaller cluster in the pair
		 */
		IReactant* smaller;

		/**
		 * The fraction rate at this size
		 */
		double fractionRate;

		/**
		 * All the coefficient needed to compute the RHS methods.
		 *
		 * 0-1 is for Xe_1
		 * 2-5 is for larger
		 * 6-9 is for smaller
		 */
		Array<double, 10> coefs;

		//! The constructor
		ReSolutionBase(IReactant* _larger, IReactant* _smaller, double _rate,
				Array<double, 10> _coefs) :
				larger(_larger), smaller(_smaller), fractionRate(_rate) {
			for (int i = 0; i < 10; i++) {
				coefs[i] = _coefs[i];
			}
		}

		/**
		 * Default constructor disallowed.
		 */
		ReSolutionBase() = delete;
		ReSolutionBase(const ReSolutionBase& other) :
				larger(other.larger), smaller(other.smaller), fractionRate(
						other.fractionRate) {
			for (int i = 0; i < 10; i++) {
				coefs[i] = other.coefs[i];
			}
		}
	};

	//! The map containing the information of who re-solves into what.
	std::vector<ReSolutionBase> sizeVec;

	//! The trap-mutation rate
	double resolutionRate;

	//! The minimum size at which the re-solution starts
	int minSize;

public:

	/**
	 * The constructor
	 */
	ReSolutionHandler() :
			resolutionRate(0.0), minSize(0) {
	}

	/**
	 * The destructor
	 */
	~ReSolutionHandler() {
	}

	/**
	 * The initialize method has to add connectivity between the Xe clusters,
	 * It must also initialize the rates of the reactions.
	 *
	 * \see IReSolutionHandler.h
	 */
	virtual void initialize(const IReactionNetwork& network,
			double electronicStoppingPower);

	/**
	 * This method update the rate for the re-solution if the fission rate
	 * changed, it should be called when the flux changes for instance.
	 *
	 * \see IReSolutionHandler.h
	 */
	void updateReSolutionRate(double rate);

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
	virtual void computeReSolution(const IReactionNetwork& network,
			double *concOffset, double *updatedConcOffset, int xi, int xs,
			int yj = 0, int zk = 0);

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
	virtual int computePartialsForReSolution(const IReactionNetwork& network,
			double *val, int *indices, int xi, int xs, int yj = 0, int zk = 0);

	/**
	 * Get the total number of clusters in the network that can undergo re-solution.
	 *
	 * \see IReSolutionHandler.h
	 */
	virtual int getNumberOfReSoluting() const {
		return sizeVec.size();
	}

	/**
	 * Set the minimum size for a cluster to undergo re-solution.
	 *
	 * \see IReSolutionHandler.h
	 */
	virtual void setMinSize(int size) {
		minSize = size;
	}

	/**
	 * Get the minimum size for a cluster to undergo re-solution.
	 *
	 * \see IReSolutionHandler.h
	 */
	virtual int getMinSize() const {
		return minSize;
	}

};
//end class ReSolutionHandler

} /* namespace xolotlCore */
#endif
