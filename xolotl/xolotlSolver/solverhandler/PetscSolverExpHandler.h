#ifndef PETSCSOLVEREXPHANDLER_H
#define PETSCSOLVEREXPHANDLER_H

// Includes
#include "PetscSolverHandler.h"
#include <Kokkos_Core.hpp>
#include <experimental/NEReactionNetwork.h>
#include <neclusters/NEClusterReactionNetwork.h>
#include <xolotlPerf/dummy/DummyHandlerRegistry.h>

namespace xolotlSolver {

/**
 * This class is a subclass of PetscSolverHandler and implement all the methods needed
 * to solve the ADR equations in 0D for the experimental network,
 * using PETSc from Argonne National Laboratory.
 */
class PetscSolverExpHandler: public PetscSolverHandler {
public:
	using NetworkType =
        xolotlCore::experimental::NEReactionNetwork;

    using ConcentrationsView = typename NetworkType::ConcentrationsView;
	using FluxesView = typename NetworkType::FluxesView;
	using SparseFillMap = typename NetworkType::SparseFillMap;
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

	// The xenon id for flux
	int xeId;

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
					*(std::make_shared<xolotlCore::NEClusterReactionNetwork>(
							std::make_shared<xolotlPerf::DummyHandlerRegistry>()))), expNetwork(
					_network), xeId(-1) {
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
	xolotlCore::experimental::NEReactionNetwork& getExpNetwork() const
			override {
		return expNetwork;
	}

};
//end class PetscSolverExpHandler

} /* end namespace xolotlSolver */

#include <solverhandler/PetscSolverExpHandler.inl>

#endif
