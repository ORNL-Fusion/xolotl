#ifndef ISORETDIFFUSIONHANDLER_H
#define ISORETDIFFUSIONHANDLER_H

// Includes
#include <memory>

#include <xolotl/core/network/IReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace modified
{
/**
 * Realizations of this interface are responsible for all the physical parts
 * for the anisotropic diffusion of mobile cluster due to temperature gradients.
 * The solver call these methods to handle the Soret diffusion.
 *
 * Note: We assume the isotropic diffusion will be ON if this is used so no
 * need to initialize ofill.
 */
class ISoretDiffusionHandler
{
public:
	/**
	 * The destructor
	 */
	virtual ~ISoretDiffusionHandler()
	{
	}

	/**
	 * Initialize which clusters are diffusing.
	 *
	 * @param network The network
	 * @param ofill Connectivity map for the off-diagonal part of the Jacobian
	 * @param dfill Connectivity map for the diagonal part of the Jacobian
	 * @param grid The X grid
	 */
	virtual void
	initialize(const network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& ofill,
		network::IReactionNetwork::SparseFillMap& dfill,
		std::vector<double> grid) = 0;

	/**
	 * This operation sets the heat flux at the surface.
	 *
	 * @param flux The heat flux at the surface
	 */
	virtual void
	setHeatFlux(double flux) = 0;

	/**
	 * This operation sets the heat conductivity to use in the equation.
	 *
	 * @param coef The heat conductivity
	 */
	virtual void
	setHeatConductivity(double cond) = 0;

	/**
	 * This operation sets the surface position.
	 *
	 * @param surfacePos The surface location
	 */
	virtual void
	updateSurfacePosition(int surfacePos) = 0;

	/**
	 * This operation sets the interface location.
	 *
	 * @param loc The interface location
	 */
	virtual void
	setInterfaceLocation(double loc) = 0;

	/**
	 * Compute the flux due to the diffusion for all the clusters that are
	 * diffusing, given the space parameters. This method is called by the
	 * RHSFunction from the solver.
	 *
	 * @param network The network
	 * @param concVector The pointer to the pointer of arrays of concentration
	 * at middle/ left/right/bottom/top/front/back grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at
	 * the grid point where the diffusion is computed
	 * @param hxLeft The step size on the left side of the point in the x
	 * direction
	 * @param hxRight The step size on the right side of the point in the x
	 * direction
	 * @param ix The position on the x grid
	 * @param sy The space parameter, depending on the grid step size in the y
	 * direction
	 * @param iy The position on the y grid
	 * @param sz The space parameter, depending on the grid step size in the z
	 * direction
	 * @param iz The position on the z grid
	 */
	virtual void
	computeDiffusion(network::IReactionNetwork& network, double** concVector,
		double* updatedConcOffset, double hxLeft, double hxRight, int ix,
		double sy = 0.0, int iy = 0, double sz = 0.0, int iz = 0) const = 0;

	/**
	 * Compute the partials due to the diffusion of all the diffusing clusters
	 * given the space parameters. This method is called by the RHSJacobian from
	 * the solver.
	 *
	 * @param network The network
	 * @param concVector The pointer to the pointer of arrays of concentration
	 * at middle/ left/right/bottom/top/front/back grid points
	 * @param val The pointer to the array that will contain the values of
	 * partials for the diffusion
	 * @param indices The pointer to the array that will contain the indices of
	 * the diffusing clusters in the network
	 * @param hxLeft The step size on the left side of the point in the x
	 * direction
	 * @param hxRight The step size on the right side of the point in the x
	 * direction
	 * @param ix The position on the x grid
	 * @param sy The space parameter, depending on the grid step size in the y
	 * direction
	 * @param iy The position on the y grid
	 * @param sz The space parameter, depending on the grid step size in the z
	 * direction
	 * @param iz The position on the z grid
	 *
	 * @return True if clusters diffused
	 */
	virtual bool
	computePartialsForDiffusion(network::IReactionNetwork& network,
		double** concVector, double* val, int* indices, double hxLeft,
		double hxRight, int ix, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) const = 0;
};
// end class ISoretDiffusionHandler

} /* end namespace modified */
} /* end namespace core */
} /* end namespace xolotl */
#endif
