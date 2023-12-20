#pragma once

// Includes
#include <xolotl/core/diffusion/DiffusionHandler.h>

namespace xolotl
{
namespace core
{
namespace diffusion
{
/**
 * This class is a subclass of DiffusionHandler for the isotropic diffusion of
 * mobile clusters in 1D.
 */
class Diffusion1DHandler : public DiffusionHandler
{
private:
	//! The vector to know which clusters are diffusing where
	std::vector<std::vector<bool>> diffusionGrid;

	//! Device copy of diffusion grid
	Kokkos::View<int**> diffusGrid;

protected:
	void
	syncDiffusionGrid();

public:
	//! The Constructor
	Diffusion1DHandler(double threshold) : DiffusionHandler(threshold)
	{
	}

	//! The Destructor
	~Diffusion1DHandler()
	{
	}

	/**
	 * \see IDiffusionHandler.h
	 */
	void
	initializeDiffusionGrid(
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny = 0, double hy = 0.0,
		int ys = 0, int nz = 0, double hz = 0.0, int zs = 0) override;

	/**
	 * Compute the flux due to the diffusion for all the cluster that are
	 * diffusing, given the space parameters. This method is called by the
	 * RHSFunction from the solver.
	 *
	 * If D is the diffusion coefficient, and C_l, C_r, C_m the left, right, and
	 * middle concentration of this cluster, the value of the flux is:
	 *
	 * D * (2.0 / [a * (a + b)]) * (C_l + [a/b] * C_r - [1.0 + (a/b)] * C_m)
	 *
	 * \see IDiffusionHandler.h
	 */
	void
	computeDiffusion(network::IReactionNetwork& network,
		const StencilConcArray& concVector,
		Kokkos::View<double*> updatedConcOffset, double hxLeft, double hxRight,
		int ix, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) const override;

	/**
	 * Compute the partials due to the diffusion of all the diffusing clusters
	 * given the space parameters. This method is called by the RHSJacobian from
	 * the solver.
	 *
	 * Using the same notation as for computeDiffusion, the partial derivative
	 * on the left grid point should be:
	 *
	 * D * (2.0 / [a * (a + b)])
	 *
	 * on the right grid point:
	 *
	 * D * (2.0 / [b * (a + b)])
	 *
	 * and on the middle grid point:
	 *
	 * - D * (2.0 / [a*b])
	 *
	 * \see IDiffusionHandler.h
	 */
	void
	computePartialsForDiffusion(network::IReactionNetwork& network,
		Kokkos::View<double*> val, double hxLeft, double hxRight, int ix,
		double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) const override;
};
// end class Diffusion1DHandler

} /* end namespace diffusion */
} /* end namespace core */
} /* end namespace xolotl */
