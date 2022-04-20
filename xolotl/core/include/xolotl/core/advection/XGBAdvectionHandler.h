#ifndef XGBADVECTIONHANDLER_H
#define XGBADVECTIONHANDLER_H

// Includes
#include <xolotl/core/advection/AdvectionHandler.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace advection
{
/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 * It represents the advection (drift) toward grain boundaries (GB) in the
 * X direction, parallel to the surface.
 */
class XGBAdvectionHandler : public AdvectionHandler
{
public:
	//! The Constructor
	XGBAdvectionHandler() : AdvectionHandler()
	{
	}

	//! The Destructor
	~XGBAdvectionHandler()
	{
	}

	/**
	 * This function initialize the list of clusters that will move through
	 * advection for grain boundaries.
	 *
	 * \see IAdvectionHandler.h
	 */
	void
	initialize(network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& ofillMap) override;

	/**
	 * The surface advection handler is in charge of initializing the grid for
	 * all the handlers, doesn't do anything here.
	 *
	 * \see IAdvectionHandler.h
	 */
	void
	initializeAdvectionGrid(std::vector<IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny = 1, double hy = 0.0,
		int ys = 0, int nz = 1, double hz = 0.0, int zs = 0) override
	{
		return;
	}

	/**
	 * \see IAdvectionHandler.h
	 */
	void
	computeAdvection(network::IReactionNetwork& network,
		const plsm::SpaceVector<double, 3>& pos, double** concVector,
		double* updatedConcOffset, double hxLeft, double hxRight, int ix,
		double hy = 0.0, int iy = 0, double hz = 0.0,
		int iz = 0) const override;

	void
	computeAdvection(network::IReactionNetwork& network,
		const plsm::SpaceVector<double, 3>& pos,
		const StencilConcArray& concVector,
		Kokkos::View<double*> updatedConcOffset, double hxLeft, double hxRight,
		int ix, double hy = 0.0, int iy = 0, double hz = 0.0,
		int iz = 0) const override;

	/**
	 * \see IAdvectionHandler.h
	 */
	void
	computePartialsForAdvection(network::IReactionNetwork& network, double* val,
		IdType* indices, const plsm::SpaceVector<double, 3>& pos, double hxLeft,
		double hxRight, int ix, double hy = 0.0, int iy = 0, double hz = 0.0,
		int iz = 0) const override;

	/**
	 * Compute the indices that will determine where the partial derivatives
	 * will be put in the Jacobian. This method is called by the RHSJacobian
	 * from the solver.
	 *
	 * Here we consider GB in the X direction.
	 *
	 * \see IAdvectionHandler.h
	 */
	std::array<int, 3>
	getStencilForAdvection(
		const plsm::SpaceVector<double, 3>& pos) const override;

	/**
	 * \see IAdvectionHandler.h
	 */
	bool
	isPointOnSink(const plsm::SpaceVector<double, 3>& pos) const override
	{
		// Return true if pos[0] is equal to location
		return fabs(location - pos[0]) < 0.001;
	}
};
// end class XGBAdvectionHandler

} /* end namespace advection */
} /* end namespace core */
} /* end namespace xolotl */
#endif
