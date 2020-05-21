#ifndef PETSCSOLVER2DHANDLER_H
#define PETSCSOLVER2DHANDLER_H

// Includes
#include "PetscSolverHandler.h"
#include <neclusters/NEClusterReactionNetwork.h>
#include <xolotlPerf/dummy/DummyHandlerRegistry.h>

namespace xolotlSolver {

/**
 * This class is a subclass of PetscSolverHandler and implement all the methods needed
 * to solve the ADR equations in 2D using PETSc from Argonne National Laboratory.
 */
class PetscSolver2DHandler: public PetscSolverHandler {
private:
	//! The position of the surface
	std::vector<int> surfacePosition;

public:

	/**
	 * Construct a PetscSolver2DHandler.
	 */
	PetscSolver2DHandler() = delete;

	/**
	 * Construct a PetscSolver2DHandler.
	 *
	 * @param _network The reaction network to use.
	 */
	PetscSolver2DHandler(NetworkType& _network) :
			PetscSolverHandler(
					*(std::make_shared<xolotlCore::NEClusterReactionNetwork>(
							std::make_shared<xolotlPerf::DummyHandlerRegistry>())),
					_network) {
	}

	//! The Destructor
	~PetscSolver2DHandler() {
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
	 * Compute the full Jacobian.
	 * \see ISolverHandler.h
	 */
	void computeJacobian(TS &ts, Vec &localC, Mat &J, PetscReal ftime);

	/**
	 * Get the position of the surface.
	 * \see ISolverHandler.h
	 */
	int getSurfacePosition(int j = -1, int k = -1) const {
		return surfacePosition[j];
	}

	/**
	 * Set the position of the surface.
	 * \see ISolverHandler.h
	 */
	void setSurfacePosition(int pos, int j = -1, int k = -1) {
		surfacePosition[j] = pos;

		return;
	}

};
//end class PetscSolver2DHandler

} /* end namespace xolotlSolver */
#endif
