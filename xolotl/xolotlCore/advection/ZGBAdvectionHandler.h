#ifndef ZGBADVECTIONHANDLER_H
#define ZGBADVECTIONHANDLER_H

// Includes
#include "AdvectionHandler.h"
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 * It represents the advection (drift) toward grain boundaries (GB) in the
 * Z direction, perpendicular to the surface. It has been observed that clusters
 * don't diffuse anymore once on the GB, the diffusion is thus cancelled out
 * here on the GB.
 */
class ZGBAdvectionHandler: public AdvectionHandler {
public:

	//! The Constructor
	ZGBAdvectionHandler() {}

	//! The Destructor
	~ZGBAdvectionHandler() {}

	/**
	 * This function initialize the list of clusters that will move through advection for
	 * grain boundaries.
	 *
	 * \see IAdvectionHandler.h
	 */
	void initialize(IReactionNetwork *network, int *ofill);

	/**
	 * Initialize an array of the dimension of the physical domain times the number of advecting
	 * clusters. For each location, True means the cluster is moving, False means it is not.
	 * Don't do anything here.
	 *
	 * \see IAdvectionHandler.h
	 */
	void initializeAdvectionGrid(std::vector<IAdvectionHandler *> advectionHandlers,
			std::vector<double> grid,
			int ny = 1, double hy = 0.0, int nz = 1, double hz = 0.0) {return;}

	/**
	 * Compute the flux due to the advection for all the helium clusters,
	 * given the space parameters and the position.
	 * This method is called by the RHSFunction from the PetscSolver.
	 * This method also removes the flux from diffusion of the advecting clusters
	 * on the GB.
	 *
	 * \see IAdvectionHandler.h
	 */
	void computeAdvection(IReactionNetwork *network,
			std::vector<double> &pos, double **concVector, double *updatedConcOffset,
			double hxLeft, double hxRight, int ix,
			double hy = 0.0, int iy = 0, double hz = 0.0, int iz = 0);

	/**
	 * Compute the partials due to the advection of all the helium clusters given
	 * the space parameters and the position.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 * This method also removes the partials from diffusion of the advecting clusters
	 * on the GB.
	 *
	 * \see IAdvectionHandler.h
	 */
	void computePartialsForAdvection(IReactionNetwork *network,
			double *val, int *indices, std::vector<double> &pos,
			double hxLeft, double hxRight, int ix,
			double hy = 0.0, int iy = 0, double hz = 0.0, int iz = 0);

	/**
	 * Compute the indices that will determine where the partial derivatives will
	 * be put in the Jacobian.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * Here we consider GB in the Z direction.
	 *
	 * \see IAdvectionHandler.h
	 */
	std::vector<int> getStencilForAdvection(std::vector<double> &pos);

	/**
	 * Check whether the grid point is located on the sink surface or not.
	 *
	 * \see IAdvectionHandler.h
	 */
	bool isPointOnSink(std::vector<double> &pos);

};
//end class ZGBAdvectionHandler

} /* end namespace xolotlCore */
#endif
