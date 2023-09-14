#pragma once

// Includes
#include <xolotl/factory/solver/SolverFactory.h>
#include <xolotl/solver/Solver.h>

namespace xolotl
{
namespace solver
{
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

	// For the monitors
	std::vector<std::vector<std::vector<double>>> _nSurf;
	std::vector<std::vector<std::vector<double>>> _nBulk;
	std::vector<std::vector<std::vector<double>>> _previousSurfFlux;
	std::vector<std::vector<std::vector<double>>> _previousBulkFlux;

	/**
	 * This operation configures the initial conditions of the grid in Xolotl.
	 *
	 * @param data The DM (data manager) created by PETSc
	 * @param solutionVector The solution vector that contains the PDE
	 * solution and which needs to be initialized.
	 * @param oldData The previous DM
	 * @param oldSolution The previous solution vector that contains the PDE
	 * solution.
	 */
	void
	setupInitialConditions(
		DM data, Vec solutionVector, DM oldData, Vec oldSolution);

public:
	/**
	 * Default constructor, deleted because we must construct using arguments.
	 */
	PetscSolver() = delete;

	PetscSolver(const options::IOptions& options);

	//! The Constructor
	PetscSolver(const std::shared_ptr<handler::ISolverHandler>& _solverHandler);

	//! The Destructor
	~PetscSolver();

	/**
	 * \see ISolver.h
	 */
	void
	initialize(int loop = 0, double time = 0.0, DM oldDA = nullptr,
		Vec oldC = nullptr) override;

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
