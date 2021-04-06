#ifndef DUMMYTRAPMUTATIONHANDLER_H
#define DUMMYTRAPMUTATIONHANDLER_H

#include <xolotl/core/modified/TrapMutationHandler.h>

namespace xolotl
{
namespace core
{
namespace modified
{
/**
 * This class realizes the ITrapMutationHandler interface responsible for the
 * modified trap-mutation of small helium clusters close to the surface. Here it
 * is a dummy class, meaning that it should not do anything.
 */
class DummyTrapMutationHandler : public TrapMutationHandler
{
private:
	/**
	 * \see TrapMutationHandler.h
	 */
	void
	initializeDepthSize(double temp)
	{
		// Initialize the vectors
		depthVec = {-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1};
		sizeVec = {0, 0, 0, 0, 0, 0, 0};

		// And don't do anything else
		return;
	}

public:
	/**
	 * The constructor
	 */
	DummyTrapMutationHandler()
	{
	}

	/**
	 * The Destructor
	 */
	~DummyTrapMutationHandler()
	{
	}

	/**
	 * \see ITrapMutationHandler.h
	 */
	void
	initialize(network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& dfill, IdType nx,
		IdType ny = 0, IdType nz = 0)
	{
		// Doesn't do anything
		return;
	}

	/**
	 * \see ITrapMutationHandler.h
	 */
	void
	initializeIndex1D(IdType surfacePos, network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, IdType nx, IdType xs)
	{
		// Doesn't do anything
		return;
	}

	/**
	 * \see ITrapMutationHandler.h
	 */
	void
	initializeIndex2D(std::vector<IdType> surfacePos,
		network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, IdType nx, IdType xs, IdType ny, double hy,
		IdType ys)
	{
		// Doesn't do anything
		return;
	}

	/**
	 * \see ITrapMutationHandler.h
	 */
	void
	initializeIndex3D(std::vector<std::vector<IdType>> surfacePos,
		network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, IdType nx, IdType xs, IdType ny, double hy,
		IdType ys, IdType nz, double hz, IdType zs)
	{
		// Doesn't do anything
		return;
	}

	/**
	 * \see ITrapMutationHandler.h
	 */
	void
	computeTrapMutation(network::IReactionNetwork& network, double* concOffset,
		double* updatedConcOffset, IdType xi, IdType yj = 0, IdType zk = 0)
	{
		// Doesn't do anything
		return;
	}

	/**
	 * \see ITrapMutationHandler.h
	 */
	size_t
	computePartialsForTrapMutation(network::IReactionNetwork& network,
		double* concOffset, double* val, IdType* indices, IdType xi,
		IdType yj = 0, IdType zk = 0)
	{
		// Doesn't do anything
		return 0;
	}

	/**
	 * \see ITrapMutationHandler.h
	 */
	size_t
	getNumberOfMutating() const
	{
		return 0;
	}
};
// end class DummyTrapMutationHandler

} /* namespace modified */
} /* namespace core */
} /* namespace xolotl */

#endif
