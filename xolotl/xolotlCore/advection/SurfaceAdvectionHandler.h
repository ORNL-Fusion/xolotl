#ifndef SURFACEADVECTIONHANDLER_H
#define SURFACEADVECTIONHANDLER_H

// Includes
#include "AdvectionHandler.h"
#include <experimental/PSIReactionNetwork.h>
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster, in the case
 * where cluster are drifting toward the surface.
 */
class SurfaceAdvectionHandler: public AdvectionHandler {
private:

	//! The vector to know which clusters are moving where
	std::vector<std::vector<std::vector<std::vector<bool> > > > advectionGrid;

public:

	//! The Constructor
	SurfaceAdvectionHandler() :
			AdvectionHandler() {
	}

	//! The Destructor
	~SurfaceAdvectionHandler() {
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
					override;

	/**
	 * Compute the flux due to the advection for all the helium clusters,
	 * given the space parameter hx and the position.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * If D is the diffusion coefficient, and C_r, C_m the right and middle concentration
	 * of this cluster, A the sink strength, K the Boltzmann constant, T the temperature,
	 * the value to add to the updated concentration is:
	 *
	 * [(3 * A * D) / (K * T * hx)] * [(C_r / [pos_x + hx]^4) - (C_m / (pos_x)^4)]
	 *
	 * \see IAdvectionHandler.h
	 */
	void computeAdvection(experimental::IReactionNetwork& network,
			const Point<3>& pos, double **concVector, double *updatedConcOffset,
			double hxLeft, double hxRight, int ix, double hy = 0.0, int iy = 0,
			double hz = 0.0, int iz = 0) const override;

	/**
	 * Compute the partials due to the advection of all the helium clusters given
	 * the space parameter hx and the position.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * The partial derivative on the right grid point is given by (same notation as for
	 * the computeAdvection method)
	 *
	 * (3 * A * D) / [K * T * hx * (pos_x + hx)^4]
	 *
	 * and on this grid point we have
	 *
	 * - (3 * A * D) / [K * T * hx * (pos_x)^4]
	 *
	 * \see IAdvectionHandler.h
	 */
	void computePartialsForAdvection(
			experimental::IReactionNetwork& network, double *val,
			int *indices, const Point<3>& pos, double hxLeft, double hxRight,
			int ix, double hy = 0.0, int iy = 0, double hz = 0.0,
			int iz = 0) const override;

	/**
	 * Compute the indices that will determine where the partial derivatives will
	 * be put in the Jacobian.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * For the surface advection the stencil is always the same:
	 *
	 * stencil[0] = 1 //x
	 * stencil[1] = 0 //y
	 * stencil[2] = 0 //z
	 *
	 * \see IAdvectionHandler.h
	 */
	std::array<int, 3> getStencilForAdvection(const Point<3>& pos) const
			override {
		// Always return (1, 0, 0)
		return {1, 0, 0};
	}

	/**
	 * Check whether the grid point is located on the sink surface or not.
	 *
	 * \see IAdvectionHandler.h
	 */
	bool isPointOnSink(const Point<3>& pos) const override {
		// Always return false
		return false;
	}

};
//end class SurfaceAdvectionHandler

} /* end namespace xolotlCore */
#endif
