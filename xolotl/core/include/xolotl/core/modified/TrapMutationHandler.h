#ifndef TRAPMUTATIONHANDLER_H
#define TRAPMUTATIONHANDLER_H

// Includes
#include <xolotl/core/Constants.h>
#include <xolotl/core/modified/ITrapMutationHandler.h>
#include <xolotl/core/modified/Sigma3TrapMutationHandler.h>

namespace xolotl
{
namespace core
{
namespace modified
{
/**
 * This class realizes ITrapMutationHandler interface, responsible for the
 * modified trap-mutation of small helium clusters close to the surface.
 */
class TrapMutationHandler : public ITrapMutationHandler
{
protected:
	/**
	 * This is a protected class used to store information about helium
	 * desorption.
	 */
	class Desorption
	{
	public:
		/**
		 * The size of the helium cluster that desorpts
		 */
		int size;

		/**
		 * Its Id
		 */
		std::size_t id;

		/**
		 * The portion of the reactions that is desorption and NOT trap-mutation
		 */
		double portion;

		//! The constructor
		Desorption(int s, double p) :
			size(s),
			portion(p),
			id(network::IReactionNetwork::invalidIndex())
		{
		}
	};

	/** The vector containing the different depths for the modified
	 * trap-mutation associated with the surface
	 */
	std::vector<double> depthVec;

	//! The vector containing the different vacancy size for the modified
	//! trap-mutation
	std::vector<int> sizeVec;

	//! The trap-mutation rate
	double kMutation;

	//! The disappearing rate
	double kDis;

	//! To know if we want attenuation or not
	bool attenuation;

	//! To know if we are using the reduced Jacobian method.
	bool enableReducedJacobian;

	/**
	 * The vector containing the indices of the bubbles created through modified
	 * trap-mutation for each grid point. The difference between this vector and
	 * depthVec is that this one is used for the actual computation whereas the
	 * other one is defined by the user. tmBubbles is created with the depthVec
	 * information.
	 *
	 * TODO: change name?
	 */
	using ReactantRefVector1D =
		std::vector<std::vector<std::array<std::size_t, 3>>>;
	using ReactantRefVector2D = std::vector<ReactantRefVector1D>;
	using ReactantRefVector3D = std::vector<ReactantRefVector2D>;
	ReactantRefVector3D tmBubbles;

	/**
	 * The desorption information
	 */
	Desorption desorp;

	/**
	 * This method fills two vectors to define the modified trap-mutation: for
	 * the first one, the first value corresponds to the depth at which the He1
	 * cluster undergo trap-mutation (if the value is negative it means that it
	 * doesn't TM), the second value correspond to He2, etc.; the second vector
	 * gives the size of the vacancies into which He trap-mutates. Information
	 * about desorption is also initialized here. It needs to be implemented by
	 * the daughter classes.
	 *
	 * @param temp The temperature of the system
	 */
	virtual void
	initializeDepthSize(double temp)
	{
		return;
	}

public:
	/**
	 * The constructor
	 */
	TrapMutationHandler() :
		kMutation(0.0),
		kDis(1.0),
		attenuation(true),
		desorp(0, 0.0),
		enableReducedJacobian(false)
	{
	}

	/**
	 * The destructor
	 */
	~TrapMutationHandler()
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
	virtual void
	initialize(network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& dfill, int nx, int ny = 0,
		int nz = 0);

	/**
	 * This method defines which trap-mutation is allowed at each grid point.
	 * The stored indices correspond to the HeV bubbles, and more precisely to
	 * their rank in the bubbles vector obtained with bubbles =
	 * network->getAll(heVType).
	 *
	 * \see ITrapMutationHandler.h
	 */
	virtual void
	initializeIndex1D(int surfacePos, network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs);

	/**
	 * This method defines which trap-mutation is allowed at each grid point.
	 * The stored indices correspond to the HeV bubbles, and more precisely to
	 * their rank in the bubbles vector obtained with bubbles =
	 * network->getAll(heVType).
	 *
	 * \see ITrapMutationHandler.h
	 */
	virtual void
	initializeIndex2D(std::vector<int> surfacePos,
		network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny, double hy, int ys);

	/**
	 * This method defines which trap-mutation is allowed at each grid point.
	 * The stored indices correspond to the HeV bubbles, and more precisely to
	 * their rank in the bubbles vector obtained with bubbles =
	 * network->getAll(heVType).
	 *
	 * \see ITrapMutationHandler.h
	 */
	virtual void
	initializeIndex3D(std::vector<std::vector<int>> surfacePos,
		network::IReactionNetwork& network,
		std::vector<advection::IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny, double hy, int ys,
		int nz, double hz, int zs);

	/**
	 * This method update the rate for the modified trap-mutation if the rates
	 * changed in the network, it should be called when temperature changes
	 * for instance.
	 *
	 * @param rate The largest rate in the network
	 */
	void
	updateTrapMutationRate(const double rate);

	/**
	 * This method set the boolean to remember if we want attenuation or not.
	 *
	 * @param isAttenuation True if we want attenuation
	 */
	void
	setAttenuation(bool isAttenuation);

	/**
	 * This method update the rate that makes the modified trap-mutation
	 * inefficient with time, depending on the total helium concentration.
	 *
	 * @param conc The concentration of helium
	 */
	void
	updateDisappearingRate(double conc);

	/**
	 * Compute the flux due to the modified trap-mutation for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * He_i --> (He_i)(V) + I
	 *
	 * F(He_i) = -F[(He_i)(V)] = -F(I) = -kMutation * C_(He_i)
	 *
	 * \see ITrapMutationHandler.h
	 */
	virtual void
	computeTrapMutation(network::IReactionNetwork& network, double* concOffset,
		double* updatedConcOffset, int xi, int yj = 0, int zk = 0);

	/**
	 * Compute the partials due to the modified trap-mutation for all the
	 * clusters given the position index xi. Returns the number of helium
	 * clusters that are undergoing trap-mutation at this grid point.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * He_i --> (He_i)(V) + I
	 *
	 * dF(He_i)/dC_(He_i) = -dF[(He_i)(V)]/dC_(He_i) = -dF(I)/dC_(He_i)
	 * 		= -kMutation
	 *
	 * \see ITrapMutationHandler.h
	 */
	virtual int
	computePartialsForTrapMutation(network::IReactionNetwork& network,
		double* concOffset, double* val, int* indices, int xi, int yj = 0,
		int zk = 0);

	/**
	 * Get the total number of clusters in the network that can undergo trap
	 * mutation.
	 *
	 * @return The number of clusters
	 */
	virtual int
	getNumberOfMutating() const
	{
		return sizeVec.size();
	}

	/**
	 * Set the enable reduced Jacobian method.
	 *
	 * @param The boolean
	 */
	void
	setEnableReduceJacobian(bool reduced)
	{
		enableReducedJacobian = reduced;
		return;
	}
};
// end class TrapMutationHandler

} /* namespace modified */
} /* namespace core */
} /* namespace xolotl */
#endif
