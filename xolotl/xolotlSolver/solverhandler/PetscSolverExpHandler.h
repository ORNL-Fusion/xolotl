#ifndef PETSCSOLVEREXPHANDLER_H
#define PETSCSOLVEREXPHANDLER_H

// Includes
#include "PetscSolverHandler.h"
#include <Kokkos_Core.hpp>
#include <experimental/PSIReactionNetwork.h>
#include <psiclusters/PSIClusterReactionNetwork.h>
#include <xolotlPerf/dummy/DummyHandlerRegistry.h>

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
	xolotlCore::experimental::PSIReactionNetwork<TImpl>;

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

	/**
	 * The current temperature on the grid.
	 */
	std::vector<double> temperature;

	// The helium id for flux
	int heId;
	int dId;
	int tId;

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
			PetscSolverHandler(
					*(std::make_shared<xolotlCore::PSIClusterReactionNetwork>(
							std::make_shared<xolotlPerf::DummyHandlerRegistry>()))), expNetwork(
					_network), heId(-1), dId(-1), tId(-1) {
	}

	//! The Destructor
	~PetscSolverExpHandler() {
	}

	/**
	 * Create everything needed before starting to solve.
	 * \see ISolverHandler.h
	 */
	void createSolverContext(DM &da) override;

	/**
	 * Initialize the concentration solution vector.
	 * \see ISolverHandler.h
	 */
	void initializeConcentration(DM &da, Vec &C) override;

	/**
	 * Compute the new concentrations for the RHS function given an initial
	 * vector of concentrations. Apply the diffusion, advection and all the reactions.
	 * \see ISolverHandler.h
	 */
	void updateConcentration(TS &ts, Vec &localC, Vec &F, PetscReal ftime)
			override;

	/**
	 * Compute the off-diagonal part of the Jacobian which is related to cluster's motion.
	 * \see ISolverHandler.h
	 */
	void computeOffDiagonalJacobian(TS &ts, Vec &localC, Mat &J,
			PetscReal ftime) override;

	/**
	 * Compute the diagonal part of the Jacobian which is related to cluster reactions.
	 * \see ISolverHandler.h
	 */
	void computeDiagonalJacobian(TS &ts, Vec &localC, Mat &J, PetscReal ftime)
			override;

	/**
	 * Get the position of the surface.
	 * \see ISolverHandler.h
	 */
	int getSurfacePosition(int j = -1, int k = -1) const override {
		return 0;
	}

	/**
	 * Set the position of the surface.
	 * \see ISolverHandler.h
	 */
	void setSurfacePosition(int pos, int j = -1, int k = -1) override {
		return;
	}

	/**
	 * Get the network.
	 * \see ISolverHandler.h
	 */
	xolotlCore::experimental::PSIReactionNetwork<
			xolotlCore::experimental::PSIFullSpeciesList>& getExpNetwork() const
			override {
		return expNetwork;
	}

};
//end class PetscSolverExpHandler

} /* end namespace xolotlSolver */

#include <solverhandler/PetscSolverExpHandler.inl>

#endif
