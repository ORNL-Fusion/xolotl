#ifndef DESORPTIONHANDLER_H
#define DESORPTIONHANDLER_H

// Includes
#include <IDesorptionHandler.h>
#include <Constants.h>

namespace xolotlCore {

/**
 * This class realizes IDesorptionHandler interface, responsible for the desorption
 *  of hydrogen clusters at the surface.
 */
class DesorptionHandler: public IDesorptionHandler {
protected:

	//! The number of clusters that can desorb
	int nCluster;

	//! The recombination rate
	double kRecombination;

	//! The equilibrium concentration
	double equilibriumConc;

	/**
	 * The vector containing the indices of the cluster undergoing desorption for each grid point.
	 */
	using ReactantRefVector1D = std::vector<IReactant::RefVector>;
	using ReactantRefVector2D = std::vector<ReactantRefVector1D>;
	using ReactantRefVector3D = std::vector<ReactantRefVector2D>;
	ReactantRefVector3D hClusters;

public:

	/**
	 * The constructor
	 */
	DesorptionHandler() :
			nCluster(2), kRecombination(0.5), equilibriumConc(1.0) {
	}

	/**
	 * The destructor
	 */
	~DesorptionHandler() {
	}

	/**
	 * The initialize method generates an empty vector for desorption.
	 *
	 * \see IDesorptionHandler.h
	 */
	virtual void initialize(std::vector<double> grid, int ny = 0, double hy = 0.0,
			int nz = 0, double hz = 0.0);

	/**
	 * This method defines which cluster desorbs at each grid point.
	 *
	 * \see IDesorptionHandler.h
	 */
	virtual void initializeIndex1D(int surfacePos, const IReactionNetwork& network,
			std::vector<double> grid);

	/**
	 * This method defines which cluster desorbs at each grid point.
	 *
	 * \see IDesorptionHandler.h
	 */
	virtual void initializeIndex2D(std::vector<int> surfacePos,
			const IReactionNetwork& network, std::vector<double> grid, int ny,
			double hy);

	/**
	 * This method defines which cluster desorbs at each grid point.
	 *
	 * \see IDesorptionHandler.h
	 */
	virtual void initializeIndex3D(std::vector<std::vector<int> > surfacePos,
			const IReactionNetwork& network, std::vector<double> grid, int ny,
			double hy, int nz, double hz);

	/**
	 * Compute the flux due to the desorption for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * H_i --> 0
	 *
	 * F(H_i) = -kRecombination * (C_(H_i))^2
	 *
	 * if C_(H_i) > C_eq
	 *
	 * \see IDesorptionHandler.h
	 */
	void computeDesorption(double *concOffset, double *updatedConcOffset,
			int xi, int xs, int yj = 0, int zk = 0);

	/**
	 * Compute the partials due to the desorption for all the
	 * clusters given the position index xi. Returns the number of
	 * clusters that are undergoing desorption at this grid point.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * H_i --> 0
	 *
	 * dF(H_i)/dC_(H_i) = - 2 * kRecombination * C_(H_i)
	 *
	 * if C_(H_i) > C_eq
	 *
	 * \see IDesorptionHandler.h
	 */
	int computePartialsForDesorption(double *val, int *indices, int xi, int xs,
			int yj = 0, int zk = 0);

	/**
	 * Get the total number of clusters in the network that can undergo desorption.
	 *
	 * @return The number of clusters
	 */
	virtual int getNumberOfDesorbing() const {
		return nCluster;
	}

};
//end class DesorptionHandler

} /* namespace xolotlCore */
#endif
