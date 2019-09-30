#ifndef DUMMYADVECTIONHANDLER_H
#define DUMMYADVECTIONHANDLER_H

// Includes
#include "SurfaceAdvectionHandler.h"

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile cluster. Here it is a dummy
 * class which means that it is not actually doing anything.
 */
class DummyAdvectionHandler: public SurfaceAdvectionHandler {

public:

	//! The Constructor
	DummyAdvectionHandler() :
			SurfaceAdvectionHandler() {
	}

	//! The Destructor
	~DummyAdvectionHandler() {
	}

	/**
	 * This function initialize the list of clusters that will move through advection. For the
	 * dummy class we don't want any cluster to advect, so this class only clears the vector
	 * and doesn't fill them.
	 *
	 * @param network The network
	 * @param ofillMap Map of connectivity for advecting clusters.
	 */
	void initialize(const IReactionNetwork& network,
			IReactionNetwork::SparseFillMap& ofillMap) override {
		// Clear the index and sink strength vectors
		advectingClusters.clear();
		sinkStrengthVector.clear();

		// Return now to leave them empty
		return;
	}

	/**
	 * Initialize an array of the dimension of the physical domain times the number of advecting
	 * clusters. For each location, True means the cluster is moving, False means it is not.
	 *
	 * \see IAdvectionHandler.h
	 */
	void initializeAdvectionGrid(
			std::vector<IAdvectionHandler *> advectionHandlers,
			std::vector<double> grid, int nx, int xs, int ny = 1, double hy =
					0.0, int ys = 0, int nz = 1, double hz = 0.0, int zs = 0)
					override {
		// Doesn't do anything
		return;
	}

	/**
	 * Compute the flux due to the advection for all the helium clusters,
	 * given the space parameter hx and the position.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * \see IAdvectionHandler.h
	 */
	void computeAdvection(const IReactionNetwork& network, const Point<3>& pos,
			double **concVector, double *updatedConcOffset, double hxLeft,
			double hxRight, int ix, double hy = 0.0, int iy = 0,
			double hz = 0.0, int iz = 0) const override {
		// Doesn't do anything
		return;
	}

	/**
	 * Compute the partials due to the advection of all the helium clusters given
	 * the space parameter hx and the position.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * \see IAdvectionHandler.h
	 */
	void computePartialsForAdvection(const IReactionNetwork& network,
			double *val, int *indices, const Point<3>& pos, double hxLeft,
			double hxRight, int ix, double hy = 0.0, int iy = 0,
			double hz = 0.0, int iz = 0) const override {
		// Doesn't do anything
		return;
	}

};
//end class DummyAdvectionHandler

} /* end namespace xolotlCore */
#endif
