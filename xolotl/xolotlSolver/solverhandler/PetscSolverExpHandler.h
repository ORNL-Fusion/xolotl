#ifndef PETSCSOLVEREXPHANDLER_H
#define PETSCSOLVEREXPHANDLER_H

// Includes
#include "PetscSolverHandler.h"
#include <Kokkos_Core.hpp>
#include <experimental/ReactionNetwork.h>

namespace xolotlSolver {

/**
 * This class is a subclass of PetscSolverHandler and implement all the methods needed
 * to solve the ADR equations in 0D for the experimental network,
 * using PETSc from Argonne National Laboratory.
 */
template<typename TImpl>
class PetscSolverExpHandler: public PetscSolverHandler {
public:
	using NetworkType =
	xolotlCore::experimental::ReactionNetwork<TImpl>;

	using ConcentrationsView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
	using FluxesView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
	using SparseFillMap = std::unordered_map<int, std::vector<int>>;
protected:

	//! The original network created from the network loader.
	NetworkType& expNetwork;

	//! Partial derivatives for all reactions at one grid point.
	Kokkos::View<double*> expVals;

	//! Map of connectivities
	SparseFillMap dfill;

public:

	/**
	 * Construct a PetscSolverExpHandler.
	 *
	 * @param _network The reaction network to use.
	 */
	PetscSolverExpHandler() = delete;

	/**
	 * Construct a PetscSolverExpHandler.
	 *
	 * @param _network The reaction network to use.
	 */
	PetscSolverExpHandler(NetworkType& _network) :
			expNetwork(_network) {
	}

	//! The Destructor
	~PetscSolverExpHandler() {
	}

	/**
	 * Create everything needed before starting to solve.
	 * \see ISolverHandler.h
	 */
	void createSolverContext(DM &da);

	/**
	 * Initialize the concentration solution vector.
	 * \see ISolverHandler.h
	 */
	void initializeConcentration(DM &da, Vec &C);

	/**
	 * Compute the new concentrations for the RHS function given an initial
	 * vector of concentrations. Apply the diffusion, advection and all the reactions.
	 * \see ISolverHandler.h
	 */
	void updateConcentration(TS &ts, Vec &localC, Vec &F, PetscReal ftime);

	/**
	 * Compute the off-diagonal part of the Jacobian which is related to cluster's motion.
	 * \see ISolverHandler.h
	 */
	void computeOffDiagonalJacobian(TS &ts, Vec &localC, Mat &J,
			PetscReal ftime);

	/**
	 * Compute the diagonal part of the Jacobian which is related to cluster reactions.
	 * \see ISolverHandler.h
	 */
	void computeDiagonalJacobian(TS &ts, Vec &localC, Mat &J, PetscReal ftime);

};
//end class PetscSolverExpHandler

} /* end namespace xolotlSolver */
#endif
