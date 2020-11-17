#ifndef DIFFUSION2DHANDLER_H
#define DIFFUSION2DHANDLER_H

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
 * mobile clusters in 2D.
 */
class Diffusion2DHandler : public DiffusionHandler
{
private:
	//! The vector to know which clusters are diffusing where
	std::vector<std::vector<std::vector<bool>>> diffusionGrid;

public:
	//! The Constructor
	Diffusion2DHandler(double threshold) : DiffusionHandler(threshold)
	{
	}

	//! The Destructor
	~Diffusion2DHandler()
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
	 * If D is the diffusion coefficient, and C_l, C_r, C_b, C_t, C_m the left,
	 * right, bottom, top, and middle concentration of this cluster, the value
	 * of the flux is:
	 *
	 * D * [(2.0 / [a * (a + b)]) * (C_l + [a/b] * C_r - [1.0 + (a/b)] * C_m)
	 *     + sy * (C_b + C_t - 2*C_m)]
	 *
	 * \see IDiffusionHandler.h
	 */
	void
	computeDiffusion(network::IReactionNetwork& network, double** concVector,
		double* updatedConcOffset, double hxLeft, double hxRight, int ix,
		double sy = 0.0, int iy = 0, double sz = 0.0,
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
	 * on the bottom and top grid points:
	 *
	 * D * sy
	 *
	 * and on the middle grid point:
	 *
	 * - 2.0 * D * [(1.0 / [a*b]) + sy]
	 *
	 * On the directly neighboring grid points the partial will be D * si with D
	 * the diffusion coefficient and the si the associated space parameter, and
	 * on this grid point the value of the partial will be -2 * D * (sx + sy).
	 *
	 * \see IDiffusionHandler.h
	 */
	void
	computePartialsForDiffusion(network::IReactionNetwork& network, double* val,
		int* indices, double hxLeft, double hxRight, int ix, double sy = 0.0,
		int iy = 0, double sz = 0.0, int iz = 0) const override;
};
// end class Diffusion2DHandler

} /* end namespace diffusion */
} /* end namespace core */
} /* end namespace xolotl */
#endif
