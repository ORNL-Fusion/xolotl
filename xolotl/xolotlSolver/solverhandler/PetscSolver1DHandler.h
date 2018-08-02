#ifndef PETSCSOLVER1DHANDLER_H
#define PETSCSOLVER1DHANDLER_H

// Includes
#include "PetscSolverHandler.h"

namespace xolotlSolver {

/**
 * This class is a subclass of PetscSolverHandler and implement all the methods needed
 * to solve the ADR equations in 1D using PETSc from Argonne National Laboratory.
 */
class PetscSolver1DHandler: public PetscSolverHandler {
private:
	//! The position of the surface
	int surfacePosition;

	/**
	 * The last temperature on the grid. It is a vector to keep the temperature at each
	 * grid point but we know the temperature changes with depth only.
	 */
	std::vector<double> lastTemperature;

public:

	/**
	 * Construct a PetscSolver1DHandler.
	 */
	PetscSolver1DHandler() = delete;

	/**
	 * Construct a PetscSolver1DHandler.
	 *
	 * @param _network The reaction network to use.
	 */
	PetscSolver1DHandler(xolotlCore::IReactionNetwork& _network) :
			PetscSolverHandler(_network), surfacePosition(0) {
	}

	//! The Destructor
	~PetscSolver1DHandler() {
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

	/**
	 * Get the position of the surface.
	 * \see ISolverHandler.h
	 */
	int getSurfacePosition(int j = -1, int k = -1) const {
		return surfacePosition;
	}

	/**
	 * Set the position of the surface.
	 * \see ISolverHandler.h
	 */
	void setSurfacePosition(int pos, int j = -1, int k = -1) {
		// Update the network rates
		network.addGridPoints(pos - surfacePosition);
		// Update the size of the last temperature vector
		if (pos > surfacePosition) {
			// Removing grid points
			while (pos > surfacePosition) {
				lastTemperature.erase(lastTemperature.begin(),
						lastTemperature.begin() + 1);
				surfacePosition++;
			}
		} else {
			// Add grid point
			while (pos < surfacePosition) {
				lastTemperature.emplace(lastTemperature.begin(), 0.0);
				surfacePosition--;
			}
		}
	}

};
//end class PetscSolver1DHandler

} /* end namespace xolotlSolver */
#endif
