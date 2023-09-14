#ifndef DUMMYDIFFUSIONHANDLER_H
#define DUMMYDIFFUSIONHANDLER_H

// Includes
#include <xolotl/core/diffusion/DiffusionHandler.h>

namespace xolotl
{
namespace core
{
namespace diffusion
{
/**
 * This class realizes the IDiffusionHandler interface responsible for all
 * the physical parts for the diffusion of mobile cluster. Here it is a dummy
 * class, meaning that it should not do anything.
 */
class DummyDiffusionHandler : public DiffusionHandler
{
public:
	//! The Constructor
	DummyDiffusionHandler(double threshold) : DiffusionHandler(threshold)
	{
	}

	//! The Destructor
	~DummyDiffusionHandler()
	{
	}

	/**
	 * Initialize the off-diagonal part of the Jacobian. If this step is skipped
	 * it won't be possible to set the partials for the diffusion.
	 *
	 * We don't want any cluster to diffuse, so nothing is set to 1 in ofill,
	 * and no index is added to the vector.
	 *
	 * \see IDiffusionHandler.h
	 */
	void
	initialize(network::IReactionNetwork&, std::vector<RowColPair>&) override
	{
		// Clear the index vector
		diffusingClusters.clear();

		// And don't do anything else
		return;
	}

	/**
	 * Initialize an array of the dimension of the physical domain times the
	 * number of diffusion clusters.
	 *
	 * \see IDiffusionHandler.h
	 */
	void
	initializeDiffusionGrid(
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny = 0, double hy = 0.0,
		int ys = 0, int nz = 0, double hz = 0.0, int zs = 0) override
	{
		// Don't do anything
		return;
	}

	/**
	 * Compute the flux due to the diffusion for all the cluster that are
	 * diffusing, given the space parameters. This method is called by the
	 * RHSFunction from the PetscSolver.
	 *
	 * Here it won't do anything because it is a dummy class.
	 *
	 * \see IDiffusionHandler.h
	 */
	void
	computeDiffusion(network::IReactionNetwork& network,
		const StencilConcArray& concVector,
		Kokkos::View<double*> updatedConcOffset, double hxLeft, double hxRight,
		int ix, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) const override
	{
		return;
	}

	/**
	 * Compute the partials due to the diffusion of all the diffusing clusters
	 * given the space parameters. This method is called by the RHSJacobian from
	 * the PetscSolver.
	 *
	 * Here it won't do anything because it is a dummy class.
	 *
	 * \see IDiffusionHandler.h
	 */
	void
	computePartialsForDiffusion(network::IReactionNetwork& network, double* val,
		IdType* indices, double hxLeft, double hxRight, int ix, double sy = 0.0,
		int iy = 0, double sz = 0.0, int iz = 0) const override
	{
		return;
	}
};
// end class DummyDiffusionHandler

} /* end namespace diffusion */
} /* end namespace core */
} /* end namespace xolotl */
#endif
