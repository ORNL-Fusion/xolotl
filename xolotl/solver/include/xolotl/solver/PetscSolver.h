#pragma once

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
 * This class realizes the Solver interface to solve the
 * diffusion-reaction problem with the PETSc solvers from Argonne
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
	 * Timer for rhsFunction
	 */
	std::shared_ptr<perf::ITimer> rhsFunctionTimer;

	/**
	 * Timer for rhsJacobian()
	 */
	std::shared_ptr<perf::ITimer> rhsJacobianTimer;

	/**
	 * Timer for solve()
	 */
	std::shared_ptr<perf::ITimer> solveTimer;

	/**
	 * This operation configures the initial conditions of the grid in Xolotl.
	 *
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
	PetscSolver(const std::shared_ptr<handler::ISolverHandler>& _solverHandler,
		const std::shared_ptr<perf::IPerfHandler>& _perfHandler);

	//! The Destructor
	~PetscSolver();

	/**
	 * \see ISolver.h
	 */
	void
	initialize() override;

	/**
	 * \see ISolver.h
	 */
	void
	setTimes(double finalTime, double dt) override;

	/**
	 * \see ISolver.h
	 */
	std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
	getConcVector() override;

	/**
	 * \see ISolver.h
	 */
	void
	setConcVector(std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>&
			concVector) override;

	/**
	 * \see ISolver.h
	 */
	double
	getCurrentDt() override;

	/**
	 * \see ISolver.h
	 */
	void
	setCurrentTimes(double currentTime, double currentDt) override;

	/**
	 * \see ISolver.h
	 */
	void
	solve() override;

	/**
	 * \see ISolver.h
	 */
	bool
	getConvergenceStatus() override;

	/**
	 * \see ISolver.h
	 */
	void
	finalize() override;

	/**
	 * \see ISolver.h
	 */
	TS&
	getTS() override
	{
		return ts;
	}

	/**
	 * \see ISolver.h
	 */
	double
	getXolotlTime() override;

	PetscErrorCode
	rhsFunction(TS ts, PetscReal ftime, Vec C, Vec F);

	PetscErrorCode
	rhsJacobian(TS ts, PetscReal ftime, Vec C, Mat A, Mat J);
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
