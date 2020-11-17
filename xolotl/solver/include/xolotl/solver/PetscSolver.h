#ifndef PETSCSOLVER_H
#define PETSCSOLVER_H

// Includes
#include <xolotl/factory/solver/SolverFactory.h>
#include <xolotl/solver/Solver.h>

namespace xolotl
{
namespace solver
{
#ifndef CHECK_PETSC_ERROR
#define CHECK_PETSC_ERROR
/**
 * This operation checks a PETSc error code and throws an exception with given
 * error message.
 *
 * @param errorCode The PETSc error code.
 * @param errMsg The error message in the thrown exception.
 */
inline void
checkPetscError(PetscErrorCode errorCode, const char* errorMsg)
{
	if (PetscUnlikely(errorCode)) {
		throw std::runtime_error(errorMsg);
	}
}
#endif

/**
 * This class realizes the ISolver interface to solve the
 * advection-diffusion-reaction problem with the PETSc solvers from Argonne
 * National Laboratory.
 */
class PetscSolver : public Solver
{
private:
	bool petscInitializedHere{false};

	/**
	 * Distributed array that will contain the solution.
	 */
	DM da;

	/**
	 * Time stepper.
	 */
	TS ts;

	/**
	 * Global solution vector.
	 */
	Vec C;

	/**
	 * Options.
	 */
	PetscOptions petscOptions;

	/**
	 * This operation configures the initial conditions of the grid in Xolotl.
	 * @param data The DM (data manager) created by PETSc
	 * @param solutionVector The solution vector that contains the PDE
	 * solution and which needs to be initialized.
	 */
	void
	setupInitialConditions(DM data, Vec solutionVector);

public:
	/**
	 * Default constructor, deleted because we must construct using arguments.
	 */
	PetscSolver() = delete;

	PetscSolver(const options::IOptions& options);

	//! The Constructor
	PetscSolver(handler::ISolverHandler& _solverHandler,
		std::shared_ptr<perf::IPerfHandler> _perfHandler);

	//! The Destructor
	~PetscSolver();

	/**
	 * This operation sets the run-time options of the solver. The map is a set
	 * of key-value std::string pairs that are interpreted by the solver. These
	 * options may change during execution, but it is up to Solvers to monitor
	 * the map for changes and they may do so at their discretion.
	 * @param options The set of options as key-value pairs with option names
	 * for keys and associated values mapped to those keys. A relevant example
	 * is "startTime" and "0.01" where both are of type std::string.
	 */
	void
	setOptions(const std::map<std::string, std::string>& options) override;

	/**
	 * This operation sets up the mesh that will be used by the solver and
	 * initializes the data on that mesh. This operation will throw an exception
	 * of type std::runtime_error if the mesh can not be setup.
	 */
	void
	setupMesh() override;

	/**
	 * This operation performs all necessary initialization for the solver
	 * possibly including but not limited to setting up MPI and loading initial
	 * conditions. If the solver can not be initialized, this operation will
	 * throw an exception of type std::runtime_error.
	 */
	void
	initialize() override;

	/**
	 * This operation sets the wanted final time and max time step allowed.
	 *
	 * @param finalTime The wanted final time.
	 * @param dt The wanted max time step.
	 */
	void
	setTimes(double finalTime, double dt) override;

	/**
	 * This operation get the concentration vector with the ids.
	 *
	 * @return The concentration vector
	 */
	virtual std::vector<
		std::vector<std::vector<std::vector<std::pair<int, double>>>>>
	getConcVector() override;

	/**
	 * This operation sets the concentration vector in the current state of the
	 * simulation.
	 *
	 * @param The concentration vector
	 */
	virtual void
	setConcVector(std::vector<
		std::vector<std::vector<std::vector<std::pair<int, double>>>>>&
			concVector) override;

	/**
	 * Get the current dt.
	 *
	 * @return  The current time step
	 */
	virtual double
	getCurrentDt() override;

	/**
	 * Set the current time and dt.
	 *
	 * @param currentTime The time
	 * @param currentDt The current time step
	 */
	virtual void
	setCurrentTimes(double currentTime, double currentDt) override;

	/**
	 * This operation directs the Solver to perform the solve. If the solve
	 * fails, it will throw an exception of type std::runtime_error.
	 */
	void
	solve() override;

	/**
	 * This operation checks the solver convergence status to decide whether
	 * it converged or not.
	 *
	 * @return true if it converged
	 */
	bool
	getConvergenceStatus() override;

	/**
	 * This operation performs all necessary finalization for the solver
	 * including but not limited to cleaning up memory, finalizing MPI and
	 * printing diagnostic information. If the solver can not be finalized,
	 * this operation will throw an exception of type std::runtime_error.
	 */
	void
	finalize() override;

	/**
	 * This operation gets the TS.
	 */
	TS&
	getTS() override
	{
		return ts;
	}

	/**
	 * This operation returns the most recent time that solver converged.
	 */
	double
	getXolotlTime() override;
};
// end class PetscSolver
} /* namespace solver */
} /* namespace xolotl */

// Some compilers (e.g., recent versions of Intel) define __func__
// to include the namespace or class scope when compiled with the C++11
// support enabled.  Others don't.  Because PETSc's PetscFunctionBeginUser
// does a straight string comparison between what we call the function name
// and what it determines from the compiler, we need a way to provide
// either the scoped name or the unscoped name.
#if defined(__ICC) || defined(__INTEL_COMPILER)
#define Actual__FUNCT__(sname, fname) sname "::" fname
#else
#define Actual__FUNCT__(sname, fname) fname
#endif /* if it is the Intel compiler */

#endif
