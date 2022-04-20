#ifndef DUMMYADVECTIONHANDLER_H
#define DUMMYADVECTIONHANDLER_H

// Includes
#include <xolotl/core/advection/SurfaceAdvectionHandler.h>

namespace xolotl
{
namespace core
{
namespace advection
{
/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile cluster. Here it is a dummy
 * class which means that it is not actually doing anything.
 */
class DummyAdvectionHandler : public SurfaceAdvectionHandler
{
public:
	//! The Constructor
	DummyAdvectionHandler() : SurfaceAdvectionHandler()
	{
	}

	//! The Destructor
	~DummyAdvectionHandler()
	{
	}

	/**
	 * This function initialize the list of clusters that will move through
	 * advection. For the dummy class we don't want any cluster to advect, so
	 * this class only clears the vector and doesn't fill them.
	 *
	 * \see IAdvectionHandler.h
	 */
	void
	initialize(network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& ofillMap) override
	{
		// Clear the index and sink strength vectors
		advectingClusters.clear();
		sinkStrengthVector.clear();

		// Return now to leave them empty
		return;
	}

	/**
	 * \see IAdvectionHandler.h
	 */
	void
	initializeAdvectionGrid(std::vector<IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny = 1, double hy = 0.0,
		int ys = 0, int nz = 1, double hz = 0.0, int zs = 0) override
	{
		// Doesn't do anything
		return;
	}

	/**
	 * \see IAdvectionHandler.h
	 */
	void
	computeAdvection(network::IReactionNetwork& network,
		const plsm::SpaceVector<double, 3>& pos, double** concVector,
		double* updatedConcOffset, double hxLeft, double hxRight, int ix,
		double hy = 0.0, int iy = 0, double hz = 0.0, int iz = 0) const override
	{
		// Doesn't do anything
		return;
	}

	void
	computeAdvection(network::IReactionNetwork& network,
		const plsm::SpaceVector<double, 3>& pos,
		const StencilConcArray& concVector,
		Kokkos::View<double*> updatedConcOffset, double hxLeft, double hxRight,
		int ix, double hy = 0.0, int iy = 0, double hz = 0.0,
		int iz = 0) const override
	{
		return;
	}

	/**
	 * \see IAdvectionHandler.h
	 */
	void
	computePartialsForAdvection(network::IReactionNetwork& network, double* val,
		IdType* indices, const plsm::SpaceVector<double, 3>& pos, double hxLeft,
		double hxRight, int ix, double hy = 0.0, int iy = 0, double hz = 0.0,
		int iz = 0) const override
	{
		// Doesn't do anything
		return;
	}
};
// end class DummyAdvectionHandler

} /* end namespace advection */
} /* end namespace core */
} /* end namespace xolotl */
#endif
