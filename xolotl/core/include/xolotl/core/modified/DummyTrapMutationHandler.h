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
	 * Method initializing the depth vector, the size vector,
	 * and desorption information.
	 *
	 * @param temp The temperature of the system
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
	 * The initialize method has to add connectivity between the He clusters and
	 * HeV_i clusters of same number of He, and I_i, with i = 1, 2. It must also
	 * initialize the rates of the reactions and define which trap-mutation is
	 * allowed at each grid point.
	 *
	 * \see ITrapMutationHandler.h
	 */
	void
	initialize(network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& dfill, int nx, int ny = 0,
		int nz = 0)
	{
		// Doesn't do anything
		return;
	}

	/**
	 * This method defines which trap-mutation is allowed at each grid point.
	 * The stored indices correspond to the HeV bubbles, and more precisely to
	 * their rank in the bubbles vector obtained with bubbles =
	 * network->getAll(heVType).
	 *
	 * \see ITrapMutationHandler.h
	 */
	void
	initializeIndex1D(int surfacePos, network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs)
	{
		// Doesn't do anything
		return;
	}

	/**
	 * This method defines which trap-mutation is allowed at each grid point.
	 * The stored indices correspond to the HeV bubbles, and more precisely to
	 * their rank in the bubbles vector obtained with bubbles =
	 * network->getAll(heVType).
	 *
	 * \see ITrapMutationHandler.h
	 */
	void
	initializeIndex2D(std::vector<int> surfacePos,
		network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny, double hy, int ys)
	{
		// Doesn't do anything
		return;
	}

	/**
	 * This method defines which trap-mutation is allowed at each grid point.
	 * The stored indices correspond to the HeV bubbles, and more precisely to
	 * their rank in the bubbles vector obtained with bubbles =
	 * network->getAll(heVType).
	 *
	 * \see ITrapMutationHandler.h
	 */
	void
	initializeIndex3D(std::vector<std::vector<int>> surfacePos,
		network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny, double hy, int ys,
		int nz, double hz, int zs)
	{
		// Doesn't do anything
		return;
	}

	/**
	 * Compute the flux due to the modified trap-mutation for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * \see ITrapMutationHandler.h
	 */
	void
	computeTrapMutation(network::IReactionNetwork& network, double* concOffset,
		double* updatedConcOffset, int xi, int yj = 0, int zk = 0)
	{
		// Doesn't do anything
		return;
	}

	/**
	 * Compute the partials due to the modified trap-mutation for all the
	 * clusters given the position index xi. Returns the number of helium
	 * clusters that are undergoing trap-mutation at this grid point.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * \see ITrapMutationHandler.h
	 */
	int
	computePartialsForTrapMutation(network::IReactionNetwork& network,
		double* concOffset, double* val, int* indices, int xi, int yj = 0,
		int zk = 0)
	{
		// Doesn't do anything
		return 0;
	}

	/**
	 * Get the total number of clusters in the network that can undergo trap
	 * mutation.
	 *
	 * @return The number of clusters
	 */
	int
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
