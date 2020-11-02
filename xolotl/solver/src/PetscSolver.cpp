// Includes
#include <cassert>
#include <fstream>
#include <iostream>

#include <xolotl/io/XFile.h>
#include <xolotl/solver/PetscSolver.h>
#include <xolotl/solver/handler/PetscSolver0DHandler.h>
#include <xolotl/solver/handler/PetscSolver1DHandler.h>
#include <xolotl/solver/handler/PetscSolver2DHandler.h>
#include <xolotl/solver/handler/PetscSolver3DHandler.h>
#include <xolotl/util/MPIUtils.h>

/*
 C_t =  -D*C_xx + A*C_x + F(C) + R(C) + D(C) from Brian Wirth's SciDAC project.

 D*C_xx  - diffusion of He and V and I
 A*C_x   - advection of He
 F(C)    - forcing function; He being created.
 R(C)    - reaction terms   (clusters combining)
 D(C)    - dissociation terms (cluster breaking up)

 Sample Options:
 -da_grid_x <nx>						 -- number of grid points in the x
 direction -ts_max_steps <maxsteps>                -- maximum number of
 time-steps to take -ts_final_time <time>                   -- maximum time to
 compute to -ts_dt <size>							 -- initial size of the time
 step

 */

namespace xolotl
{
namespace solver
{
namespace detail
{
auto petscSolverRegistrations =
	xolotl::factory::solver::SolverFactory::RegistrationCollection<PetscSolver>(
		{"PETSc"});
}

// Timer for RHSFunction()
std::shared_ptr<perf::ITimer> RHSFunctionTimer;

// Timer for RHSJacobian()
std::shared_ptr<perf::ITimer> RHSJacobianTimer;

// Timer for solve()
std::shared_ptr<perf::ITimer> SolveTimer;

// PETSc options object
PetscOptions* petscOptions = nullptr;

// Help message
static char help[] = "Solves C_t =  -D*C_xx + A*C_x + F(C) + R(C) + D(C) from "
					 "Brian Wirth's SciDAC project.\n";

// ----- GLOBAL VARIABLES ----- //
namespace monitor
{
extern PetscErrorCode setupPetsc0DMonitor(TS);
extern PetscErrorCode setupPetsc1DMonitor(TS);
extern PetscErrorCode setupPetsc2DMonitor(TS);
extern PetscErrorCode setupPetsc3DMonitor(TS);
} // namespace monitor

void
PetscSolver::setupInitialConditions(DM da, Vec C)
{
	// Initialize the concentrations in the solution vector
	auto& solverHandler = Solver::getSolverHandler();
	solverHandler.initializeConcentration(da, C);

	return;
}

/* ------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "RHSFunction")
/*
 RHSFunction - Evaluates the right-hand-side of the nonlinear function defining
 the ODE

 Input Parameters:
 .  ts - the TS context
 .  ftime - the physical time at which the function is evaluated
 .  C - input vector
 .  ptr - optional user-defined context

 Output Parameter:
 .  F - function values
 */
/* ------------------------------------------------------------------- */
PetscErrorCode
RHSFunction(TS ts, PetscReal ftime, Vec C, Vec F, void*)
{
	// Start the RHSFunction Timer
	RHSFunctionTimer->start();

	PetscErrorCode ierr;

	// Get the local data vector from PETSc
	PetscFunctionBeginUser;
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);
	Vec localC;
	ierr = DMGetLocalVector(da, &localC);
	CHKERRQ(ierr);

	// Scatter ghost points to local vector, using the 2-step process
	// DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
	// By placing code between these two statements, computations can be
	// done while messages are in transition.
	ierr = DMGlobalToLocalBegin(da, C, INSERT_VALUES, localC);
	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, C, INSERT_VALUES, localC);
	CHKERRQ(ierr);

	// Set the initial values of F
	ierr = VecSet(F, 0.0);
	CHKERRQ(ierr);

	// Compute the new concentrations
	auto& solverHandler = Solver::getSolverHandler();
	solverHandler.updateConcentration(ts, localC, F, ftime);

	// Stop the RHSFunction Timer
	RHSFunctionTimer->stop();

	// Return the local vector
	ierr = DMRestoreLocalVector(da, &localC);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "RHSJacobian")
/*
 Compute the Jacobian entries based on IFunction() and insert them into the
 matrix
 */
PetscErrorCode
RHSJacobian(TS ts, PetscReal ftime, Vec C, Mat A, Mat J, void*)
{
	// Start the RHSJacobian timer
	RHSJacobianTimer->start();

	PetscErrorCode ierr;

	// Get the matrix from PETSc
	PetscFunctionBeginUser;
	ierr = MatZeroEntries(J);
	CHKERRQ(ierr);
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);
	Vec localC;
	ierr = DMGetLocalVector(da, &localC);
	CHKERRQ(ierr);

	// Get the complete data array
	ierr = DMGlobalToLocalBegin(da, C, INSERT_VALUES, localC);
	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, C, INSERT_VALUES, localC);
	CHKERRQ(ierr);

	// Get the solver handler
	auto& solverHandler = Solver::getSolverHandler();

	/* ----- Compute the partial derivatives for the reaction term ----- */
	solverHandler.computeJacobian(ts, localC, J, ftime);

	// Return the local vector
	ierr = DMRestoreLocalVector(da, &localC);
	CHKERRQ(ierr);

	ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);
	ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);

	if (A != J) {
		ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
		CHKERRQ(ierr);
		ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
		CHKERRQ(ierr);
	}

	//	ierr = MatView(J, PETSC_VIEWER_STDOUT_WORLD);

	// Stop the RHSJacobian timer
	RHSJacobianTimer->stop();

	PetscFunctionReturn(0);
}

PetscSolver::PetscSolver(const options::IOptions& options) :
	Solver(options,
		[&options](core::network::IReactionNetwork& network)
			-> std::shared_ptr<handler::ISolverHandler> {
			switch (options.getDimensionNumber()) {
			case 0:
				return std::make_shared<handler::PetscSolver0DHandler>(network);
			case 1:
				return std::make_shared<handler::PetscSolver1DHandler>(network);
			case 2:
				return std::make_shared<handler::PetscSolver2DHandler>(network);
			case 3:
				return std::make_shared<handler::PetscSolver3DHandler>(network);
			default:
				// The asked dimension is not good (e.g. -1, 4)
				throw std::string(
					"\nxolotlFactory: Bad dimension for the solver handler.");
			}
		})
{
	this->setCommandLineOptions(options.getPetscArg());

	RHSFunctionTimer = perfHandler->getTimer("RHSFunctionTimer");
	RHSJacobianTimer = perfHandler->getTimer("RHSJacobianTimer");
	SolveTimer = perfHandler->getTimer("SolveTimer");
}

PetscSolver::PetscSolver(handler::ISolverHandler& _solverHandler,
	std::shared_ptr<perf::IPerfHandler> _perfHandler) :
	Solver(_solverHandler, _perfHandler)
{
	RHSFunctionTimer = perfHandler->getTimer("RHSFunctionTimer");
	RHSJacobianTimer = perfHandler->getTimer("RHSJacobianTimer");
	SolveTimer = perfHandler->getTimer("SolveTimer");
}

PetscSolver::~PetscSolver()
{
}

void
PetscSolver::setOptions(const std::map<std::string, std::string>&)
{
}

void
PetscSolver::setupMesh()
{
}

inline bool
isPetscInitialized()
{
	PetscBool ret;
	PetscInitialized(&ret);
	return ret;
}

void
PetscSolver::initialize()
{
	PetscErrorCode ierr;

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Initialize program
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if (!isPetscInitialized()) {
		ierr = PetscInitialize(NULL, NULL, NULL, help);
		checkPetscError(
			ierr, "PetscSolver::initialize: PetscInitialize failed.");
		petscInitializedHere = true;
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create the solver options
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = PetscOptionsCreate(&petscOptions);
	checkPetscError(
		ierr, "PetscSolver::initialize: PetscOptionsCreate failed.");
	ierr = PetscOptionsInsertString(petscOptions, optionsString.c_str());
	checkPetscError(
		ierr, "PetscSolver::initialize: PetscOptionsInsertString failed.");
	ierr = PetscOptionsPush(petscOptions);
	checkPetscError(ierr, "PetscSolver::initialize: PetscOptionsPush failed.");

	// Check the option -snes_mf_operator
	PetscBool flagReduced;
	ierr = PetscOptionsHasName(NULL, NULL, "-snes_mf_operator", &flagReduced);
	auto mutationHandler = getSolverHandler().getMutationHandler();
	mutationHandler->setEnableReducedJacobian(flagReduced);

	// Create the solver context
	getSolverHandler().createSolverContext(da);

	/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Extract global vector from DMDA to hold solution
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = DMCreateGlobalVector(da, &C);
	checkPetscError(
		ierr, "PetscSolver::initialize: DMCreateGlobalVector failed.");

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create timestepping solver context
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	// Get the MPI communicator
	auto xolotlComm = util::getMPIComm();
	ierr = TSCreate(xolotlComm, &ts);
	checkPetscError(ierr, "PetscSolver::initialize: TSCreate failed.");
	ierr = TSSetType(ts, TSARKIMEX);
	checkPetscError(ierr, "PetscSolver::initialize: TSSetType failed.");
	ierr = TSARKIMEXSetFullyImplicit(ts, PETSC_TRUE);
	checkPetscError(
		ierr, "PetscSolver::initialize: TSARKIMEXSetFullyImplicit failed.");
	ierr = TSSetDM(ts, da);
	checkPetscError(ierr, "PetscSolver::initialize: TSSetDM failed.");
	ierr = TSSetProblemType(ts, TS_NONLINEAR);
	checkPetscError(ierr, "PetscSolver::initialize: TSSetProblemType failed.");
	ierr = TSSetRHSFunction(ts, NULL, RHSFunction, NULL);
	checkPetscError(ierr, "PetscSolver::initialize: TSSetRHSFunction failed.");
	ierr = TSSetRHSJacobian(ts, NULL, NULL, RHSJacobian, NULL);
	checkPetscError(ierr, "PetscSolver::initialize: TSSetRHSJacobian failed.");
	ierr = TSSetSolution(ts, C);
	checkPetscError(ierr, "PetscSolver::initialize: TSSetSolution failed.");

	// Read the times if the information is in the HDF5 file
	auto fileName = getSolverHandler().getNetworkName();
	double time = 0.0, deltaTime = 1.0e-12;
	if (!fileName.empty()) {
		io::XFile xfile(fileName);
		auto concGroup = xfile.getGroup<io::XFile::ConcentrationGroup>();
		if (concGroup and concGroup->hasTimesteps()) {
			auto tsGroup = concGroup->getLastTimestepGroup();
			assert(tsGroup);
			std::tie(time, deltaTime) = tsGroup->readTimes();

			// Give the values to the solver
			ierr = TSSetTime(ts, time);
			checkPetscError(ierr, "PetscSolver::initialize: TSSetTime failed.");
			ierr = TSSetTimeStep(ts, deltaTime);
			checkPetscError(
				ierr, "PetscSolver::initialize: TSSetTimeStep failed.");
		}
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set solver options
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = TSSetFromOptions(ts);
	checkPetscError(ierr, "PetscSolver::initialize: TSSetFromOptions failed.");

	// Switch on the number of dimensions to set the monitors
	int dim = getSolverHandler().getDimension();
	switch (dim) {
	case 0:
		// One dimension
		ierr = monitor::setupPetsc0DMonitor(ts);
		checkPetscError(
			ierr, "PetscSolver::initialize: setupPetsc0DMonitor failed.");
		break;
	case 1:
		// One dimension
		ierr = monitor::setupPetsc1DMonitor(ts);
		checkPetscError(
			ierr, "PetscSolver::initialize: setupPetsc1DMonitor failed.");
		break;
	case 2:
		// Two dimensions
		ierr = monitor::setupPetsc2DMonitor(ts);
		checkPetscError(
			ierr, "PetscSolver::initialize: setupPetsc2DMonitor failed.");
		break;
	case 3:
		// Three dimensions
		ierr = monitor::setupPetsc3DMonitor(ts);
		checkPetscError(
			ierr, "PetscSolver::initialize: setupPetsc3DMonitor failed.");
		break;
	default:
		throw std::string("PetscSolver Exception: Wrong number of dimensions "
						  "to set the monitors.");
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set initial conditions
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	setupInitialConditions(da, C);

	// Set the output precision for std::out
	std::cout.precision(16);

	// Pop the options
	ierr = PetscOptionsPop();
	checkPetscError(ierr, "PetscSolver::initialize: PetscOptionsPop failed.");

	return;
}

void
PetscSolver::setTimes(double finalTime, double dt)
{
	PetscErrorCode ierr;

	// Get the default values for the dt
	TSAdapt adapt;
	ierr = TSGetAdapt(ts, &adapt);
	checkPetscError(ierr, "PetscSolver::setTimes: TSGetAdapt failed.");
	PetscReal hmin, hmax;
	ierr = TSAdaptGetStepLimits(adapt, &hmin, &hmax);
	checkPetscError(
		ierr, "PetscSolver::setTimes: TSAdaptGetStepLimits failed.");
	// Set the new max value
	ierr = TSAdaptSetStepLimits(adapt, hmin, dt);
	checkPetscError(
		ierr, "PetscSolver::setTimes: TSAdaptSetStepLimits failed.");

	// Give the final time value to the solver
	ierr = TSSetMaxTime(ts, finalTime);
	checkPetscError(ierr, "PetscSolver::setTimes: TSSetMaxTime failed.");

	return;
}

std::vector<std::vector<std::vector<std::vector<std::pair<int, double>>>>>
PetscSolver::getConcVector()
{
	auto& solverHandler = Solver::getSolverHandler();

	return solverHandler.getConcVector(da, C);
}

void
PetscSolver::setConcVector(
	std::vector<std::vector<std::vector<std::vector<std::pair<int, double>>>>>&
		concVector)
{
	auto& solverHandler = Solver::getSolverHandler();

	return solverHandler.setConcVector(da, C, concVector);
}

double
PetscSolver::getCurrentDt()
{
	PetscErrorCode ierr;

	// Get the value from the solver
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	checkPetscError(ierr, "PetscSolver::getCurrentDt: TSGetTimeStep failed.");

	return currentTimeStep;
}

void
PetscSolver::setCurrentTimes(double time, double dt)
{
	PetscErrorCode ierr;

	// Give the values to the solver
	ierr = TSSetTime(ts, time);
	checkPetscError(ierr, "PetscSolver::setCurrentTimes: TSSetTime failed.");
	ierr = TSSetTimeStep(ts, dt);
	checkPetscError(
		ierr, "PetscSolver::setCurrentTimes: TSSetTimeStep failed.");

	return;
}

void
PetscSolver::solve()
{
	PetscErrorCode ierr;
	// Start the solve Timer
	SolveTimer->start();

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Solve the ODE system
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if (ts != NULL && C != NULL) {
		// Push the options for the solve
		ierr = PetscOptionsPush(petscOptions);
		checkPetscError(ierr, "PetscSolver::solve: PetscOptionsPush failed.");
		if (!petscInitializedHere) {
			// Reset the time step number
			ierr = TSSetStepNumber(ts, 0);
			checkPetscError(
				ierr, "PetscSolver::solve: Reset Step Number failed.");
			// Reset the GB location
			auto& solverHandler = Solver::getSolverHandler();
			solverHandler.initGBLocation(da, C);
		}
		// Start the PETSc Solve
		ierr = TSSolve(ts, C);
		checkPetscError(ierr, "PetscSolver::solve: TSSolve failed.");

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		 Write in a file if everything went well or not.
		 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   */
		// Check the option -check_collapse
		PetscBool flagCheck;
		ierr = PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCheck);
		checkPetscError(ierr,
			"PetscSolver::solve: PetscOptionsHasName (-check_collapse) "
			"failed.");
		if (flagCheck) {
			// Open the output file
			std::ofstream outputFile;
			outputFile.open("solverStatus.txt");

			// Get the converged reason from PETSc
			TSConvergedReason reason;
			ierr = TSGetConvergedReason(ts, &reason);
			checkPetscError(
				ierr, "PetscSolver::solve: TSGetConvergedReason failed.");

			// Write it
			if (reason == TS_CONVERGED_EVENT)
				outputFile << "collapsed" << std::endl;
			else if (reason == TS_DIVERGED_NONLINEAR_SOLVE ||
				reason == TS_DIVERGED_STEP_REJECTED)
				outputFile << "diverged" << std::endl;
			else
				outputFile << "good" << std::endl;

			outputFile.close();
		}

		// Popping the option database
		ierr = PetscOptionsPop();
		checkPetscError(ierr, "PetscSolver::solve: PetscOptionsPop failed.");
	}
	else {
		throw std::string("PetscSolver Exception: Unable to solve! Data not "
						  "configured properly.");
	}
	// Stop the timer
	SolveTimer->stop();

	return;
}

bool
PetscSolver::getConvergenceStatus()
{
	PetscErrorCode ierr;

	// Get the converged reason from PETSc
	TSConvergedReason reason;
	ierr = TSGetConvergedReason(ts, &reason);
	checkPetscError(ierr,
		"PetscSolver::getConvergenceStatus: TSGetConvergedReason failed.");

	// Check if it was sent by the user
	bool userReason = false;
	if (reason == TS_CONVERGED_USER)
		userReason = true;

	// Share the information with all the processes
	bool totalReason = false;
	auto xolotlComm = util::getMPIComm();
	MPI_Allreduce(
		&userReason, &totalReason, 1, MPI_C_BOOL, MPI_LOR, xolotlComm);

	if (totalReason)
		return false;

	return true;
}

void
PetscSolver::finalize()
{
	PetscErrorCode ierr;

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Free work space.
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = PetscOptionsDestroy(&petscOptions);
	checkPetscError(ierr, "PetscSolver::solve: PetscOptionsDestroy failed.");
	ierr = VecDestroy(&C);
	checkPetscError(ierr, "PetscSolver::solve: VecDestroy failed.");
	ierr = TSDestroy(&ts);
	checkPetscError(ierr, "PetscSolver::solve: TSDestroy failed.");
	ierr = DMDestroy(&da);
	checkPetscError(ierr, "PetscSolver::solve: DMDestroy failed.");

	if (petscInitializedHere) {
		ierr = PetscFinalize();
		checkPetscError(ierr, "PetscSolver::finalize: PetscFinalize failed.");
	}

	return;
}

double
PetscSolver::getXolotlTime()
{
	PetscErrorCode ierr;

	// The most recent time that Xolotl converged
	PetscReal CurrentXolotlTime;
	ierr = TSGetTime(ts, &CurrentXolotlTime);
	checkPetscError(ierr, "PetscSolver::getXolotlTime: TSGetTime failed.");
	return CurrentXolotlTime;
}
} /* end namespace solver */
} /* end namespace xolotl */
