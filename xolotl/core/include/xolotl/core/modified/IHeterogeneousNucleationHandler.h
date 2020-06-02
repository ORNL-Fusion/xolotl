#ifndef IHETEROGENEOUSNUCLEATIONHANDLER_H
#define IHETEROGENEOUSNUCLEATIONHANDLER_H

// Includes
#include <memory>
#include <xolotl/core/network/IReactionNetwork.h>

namespace xolotl {
namespace core {
namespace modified {

/**
 * Realizations of this interface are responsible for all the physical parts
 * for the heterogeneous nucleation of xenon clusters.
 * The solver call these methods to handle the heterogeneous nucleation.
 */
class IHeterogeneousNucleationHandler {

public:

	/**
	 * The destructor
	 */
	virtual ~IHeterogeneousNucleationHandler() {
	}

	/**
	 * The initialize method has to add connectivity between the Xe clusters,
	 * It must also initialize the rates of the reaction.
	 *
	 * @param network The network
	 * @param dfill the matrix for the reaction connectivities
	 */
	virtual void initialize(network::IReactionNetwork& network,
			network::IReactionNetwork::SparseFillMap& dfill) = 0;

	/**
	 * This method update the rate for the heterogeneous nucleation if the fission rate
	 * changed, it should be called when the flux changes for instance.
	 *
	 * @param rate The fission rate
	 */
	virtual void updateHeterogeneousNucleationRate(double rate) = 0;

	/**
	 * This method updates the fission yield.
	 *
	 * @param yield The fission yield
	 */
	virtual void setFissionYield(double yield) = 0;

	/**
	 * Compute the flux due to the heterogeneous nucleation for all the clusters,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param network The network
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the nucleation is computed
	 * @param updatedConcOffset The pointer to the array of the concentration
	 * at the grid point where the nucleation is computed used to find the
	 * next solution
	 * @param xi The index of the position on the grid in the depth direction
	 * @param xs The beginning of the grid on this process
	 * @param yj The index of the position on the grid in the Y direction
	 * @param zk The index of the position on the grid in the Z direction
	 */
	virtual void computeHeterogeneousNucleation(
			network::IReactionNetwork& network, double *concOffset,
			double *updatedConcOffset, int xi, int xs, int yj = 0,
			int zk = 0) = 0;

	/**
	 * Compute the partials due to the heterogeneous nucleation for all the
	 * clusters given the position index xi.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * @param network The network
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the nucleation is computed
	 * @param val The pointer to the array that will contain the values of
	 * partials for the nucleation
	 * @param indices The pointer to the array that will contain the indices
	 * of the clusters
	 * @param xi The index of the position on the grid in the depth direction
	 * @param xs The beginning of the grid on this process
	 * @param yj The index of the position on the grid in the Y direction
	 * @param zk The index of the position on the grid in the Z direction
	 *
	 * @return true if nucleation is happening
	 */
	virtual bool computePartialsForHeterogeneousNucleation(
			network::IReactionNetwork& network, double *concOffset,
			double *val, int *indices, int xi, int xs, int yj = 0,
			int zk = 0) = 0;

};
//end class IHeterogeneousNucleationHandler

} /* namespace modified */
} /* namespace core */
} /* namespace xolotl */
#endif
