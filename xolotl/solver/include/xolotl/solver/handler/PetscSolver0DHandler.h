#ifndef PETSCSOLVER0DHANDLER_H
#define PETSCSOLVER0DHANDLER_H

// Includes
#include <xolotl/solver/handler/PetscSolverHandler.h>

namespace xolotl
{
namespace solver
{
namespace handler
{
/**
 * This class is a subclass of PetscSolverHandler and implement all the methods
 * needed to solve the ADR equations in 0D using PETSc from Argonne National
 * Laboratory.
 */
class PetscSolver0DHandler : public PetscSolverHandler
{
public:
	/**
	 * Construct a PetscSolver0DHandler.
	 *
	 * @param _network The reaction network to use.
	 */
	PetscSolver0DHandler() = delete;

	/**
	 * Construct a PetscSolver0DHandler.
	 *
	 * @param _network The reaction network to use.
	 */
	PetscSolver0DHandler(NetworkType& _network) : PetscSolverHandler(_network)
	{
	}

	//! The Destructor
	~PetscSolver0DHandler()
	{
	}

	/**
	 * Create everything needed before starting to solve.
	 * \see ISolverHandler.h
	 */
	void
	createSolverContext(DM& da);

	/**
	 * Initialize the concentration solution vector.
	 * \see ISolverHandler.h
	 */
	void
	initializeConcentration(DM& da, Vec& C);

	/**
	 * Set the concentrations to 0.0 where the GBs are.
	 * \see ISolverHandler.h
	 */
	void
	initGBLocation(DM& da, Vec& C)
	{
		// Doesn't do anything in 0D
		return;
	}

	/**
	 * This operation get the concentration vector with the ids.
	 * \see ISolverHandler.h
	 */
	std::vector<std::vector<std::vector<std::vector<std::pair<int, double>>>>>
	getConcVector(DM& da, Vec& C);

	/**
	 * This operation sets the concentration vector in the current state of the
	 * simulation. \see ISolverHandler.h
	 */
	void
	setConcVector(DM& da, Vec& C,
		std::vector<
			std::vector<std::vector<std::vector<std::pair<int, double>>>>>&
			concVector);

	/**
	 * Compute the new concentrations for the RHS function given an initial
	 * vector of concentrations. Apply the diffusion, advection and all the
	 * reactions. \see ISolverHandler.h
	 */
	void
	updateConcentration(TS& ts, Vec& localC, Vec& F, PetscReal ftime);

	/**
	 * Compute the full Jacobian.
	 * \see ISolverHandler.h
	 */
	void
	computeJacobian(TS& ts, Vec& localC, Mat& J, PetscReal ftime);

	/**
	 * Get the position of the surface.
	 * \see ISolverHandler.h
	 */
	int
	getSurfacePosition(int j = -1, int k = -1) const
	{
		return 0;
	}

	/**
	 * Set the position of the surface.
	 * \see ISolverHandler.h
	 */
	void
	setSurfacePosition(int pos, int j = -1, int k = -1)
	{
		return;
	}
};
// end class PetscSolver0DHandler

} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */
#endif