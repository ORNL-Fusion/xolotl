// Includes
#include <cassert>
#include <PetscSolver.h>
#include <fstream>
#include <iostream>
#include <XFile.h>
#include <MPIUtils.h>

using namespace xolotlCore;

/*
 C_t =  -D*C_xx + A*C_x + F(C) + R(C) + D(C) from Brian Wirth's SciDAC project.

 D*C_xx  - diffusion of He and V and I
 A*C_x   - advection of He
 F(C)    - forcing function; He being created.
 R(C)    - reaction terms   (clusters combining)
 D(C)    - dissociation terms (cluster breaking up)

 Sample Options:
 -da_grid_x <nx>						 -- number of grid points in the x direction
 -ts_max_steps <maxsteps>                -- maximum number of time-steps to take
 -ts_final_time <time>                   -- maximum time to compute to
 -ts_dt <size>							 -- initial size of the time step

 */

namespace xolotlSolver {

//! PETSc options object
PetscOptions* petscOptions = nullptr;

//!Timer for RHSFunction()
std::shared_ptr<xolotlPerf::ITimer> RHSFunctionTimer;

//!Timer for RHSJacobian()
std::shared_ptr<xolotlPerf::ITimer> RHSJacobianTimer;

//! Help message
static char help[] =
		"Solves C_t =  -D*C_xx + A*C_x + F(C) + R(C) + D(C) from Brian Wirth's SciDAC project.\n";

// ----- GLOBAL VARIABLES ----- //
extern PetscErrorCode setupPetsc0DMonitor(TS&);
extern PetscErrorCode setupPetsc1DMonitor(TS&,
		std::shared_ptr<xolotlPerf::IHandlerRegistry>);
extern PetscErrorCode setupPetsc2DMonitor(TS&);
extern PetscErrorCode setupPetsc3DMonitor(TS&);
extern PetscErrorCode reset0DMonitor();
extern PetscErrorCode reset1DMonitor();
extern PetscErrorCode reset2DMonitor();
extern PetscErrorCode reset3DMonitor();

void PetscSolver::setupInitialConditions(DM da, Vec C) {
	// Initialize the concentrations in the solution vector
	auto& solverHandler = Solver::getSolverHandler();
	solverHandler.initializeConcentration(da, C);

	return;
}

/* ------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "RHSFunction")
/*
 RHSFunction - Evaluates the right-hand-side of the nonlinear function defining the ODE

 Input Parameters:
 .  ts - the TS context
 .  ftime - the physical time at which the function is evaluated
 .  C - input vector
 .  ptr - optional user-defined context

 Output Parameter:
 .  F - function values
 */
/* ------------------------------------------------------------------- */
PetscErrorCode RHSFunction(TS ts, PetscReal ftime, Vec C, Vec F, void *) {
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

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "RHSJacobian")
/*
 Compute the Jacobian entries based on IFunction() and insert them into the matrix
 */
PetscErrorCode RHSJacobian(TS ts, PetscReal ftime, Vec C, Mat A, Mat J,
		void *) {
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

	/* ----- Compute the off-diagonal part of the Jacobian ----- */
	solverHandler.computeOffDiagonalJacobian(ts, localC, J, ftime);

	ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);
	ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);

	/* ----- Compute the partial derivatives for the reaction term ----- */
	solverHandler.computeDiagonalJacobian(ts, localC, J, ftime);

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

PetscSolver::PetscSolver(ISolverHandler& _solverHandler,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		Solver(_solverHandler, registry) {
	RHSFunctionTimer = handlerRegistry->getTimer("RHSFunctionTimer");
	RHSJacobianTimer = handlerRegistry->getTimer("RHSJacobianTimer");
}

PetscSolver::~PetscSolver() {
}

void PetscSolver::setOptions(const std::map<std::string, std::string>&) {
}

void PetscSolver::setupMesh() {
}

void PetscSolver::initialize() {
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Initialize program
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscInitialize(&numCLIArgs, &CLIArgs, (char*) 0, help);

	return;
}

void PetscSolver::solve() {
	PetscErrorCode ierr;

	// Create the solver context
	DM da;
	getSolverHandler().createSolverContext(da);

	/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Extract global vector from DMDA to hold solution
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	Vec C;
	ierr = DMCreateGlobalVector(da, &C);
	checkPetscError(ierr, "PetscSolver::solve: DMCreateGlobalVector failed.");

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create timestepping solver context
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	// Get the MPI communicator
	auto xolotlComm = xolotlCore::MPIUtils::getMPIComm();
	TS ts;
	ierr = TSCreate(xolotlComm, &ts);
	checkPetscError(ierr, "PetscSolver::solve: TSCreate failed.");
	ierr = TSSetType(ts, TSARKIMEX);
	checkPetscError(ierr, "PetscSolver::solve: TSSetType failed.");
	ierr = TSARKIMEXSetFullyImplicit(ts, PETSC_TRUE);
	checkPetscError(ierr,
			"PetscSolver::solve: TSARKIMEXSetFullyImplicit failed.");
	ierr = TSSetDM(ts, da);
	checkPetscError(ierr, "PetscSolver::solve: TSSetDM failed.");
	ierr = TSSetProblemType(ts, TS_NONLINEAR);
	checkPetscError(ierr, "PetscSolver::solve: TSSetProblemType failed.");
	ierr = TSSetRHSFunction(ts, NULL, RHSFunction, NULL);
	checkPetscError(ierr, "PetscSolver::solve: TSSetRHSFunction failed.");
	ierr = TSSetRHSJacobian(ts, NULL, NULL, RHSJacobian, NULL);
	checkPetscError(ierr, "PetscSolver::solve: TSSetRHSJacobian failed.");
	ierr = TSSetSolution(ts, C);
	checkPetscError(ierr, "PetscSolver::solve: TSSetSolution failed.");

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set solver options
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscOptions petscOptions;
	ierr = PetscOptionsCreate(&petscOptions);
	checkPetscError(ierr, "PetscSolver::solve: PetscOptionsCreate failed.");
	ierr = PetscOptionsInsertString(petscOptions, optionsString.c_str());
	checkPetscError(ierr,
			"PetscSolver::solve: PetscOptionsInsertString failed.");
	ierr = PetscObjectSetOptions((PetscObject) ts, petscOptions);
	checkPetscError(ierr, "PetscSolver::solve: PetscObjectSetOptions failed.");
	ierr = TSSetFromOptions(ts);
	checkPetscError(ierr, "PetscSolver::solve: TSSetFromOptions failed.");

	// Read the times if the information is in the HDF5 file
	auto fileName = getSolverHandler().getNetworkName();
	double time = 0.0, deltaTime = 1.0e-12;
	if (!fileName.empty()) {

		XFile xfile(fileName);
		auto concGroup = xfile.getGroup<XFile::ConcentrationGroup>();
		if (concGroup and concGroup->hasTimesteps()) {
			auto tsGroup = concGroup->getLastTimestepGroup();
			assert(tsGroup);
			std::tie(time, deltaTime) = tsGroup->readTimes();

			// Give the values to the solver
			ierr = TSSetTime(ts, time);
			checkPetscError(ierr, "PetscSolver::solve: TSSetTime failed.");
			ierr = TSSetTimeStep(ts, deltaTime);
			checkPetscError(ierr, "PetscSolver::solve: TSSetTimeStep failed.");
		}
	}

	// Switch on the number of dimensions to set the monitors
	int dim = getSolverHandler().getDimension();
	switch (dim) {
	case 0:
		// One dimension
		ierr = setupPetsc0DMonitor(ts);
		checkPetscError(ierr,
				"PetscSolver::solve: setupPetsc0DMonitor failed.");
		break;
	case 1:
		// One dimension
		ierr = setupPetsc1DMonitor(ts, handlerRegistry);
		checkPetscError(ierr,
				"PetscSolver::solve: setupPetsc1DMonitor failed.");
		break;
	case 2:
		// Two dimensions
		ierr = setupPetsc2DMonitor(ts);
		checkPetscError(ierr,
				"PetscSolver::solve: setupPetsc2DMonitor failed.");
		break;
	case 3:
		// Three dimensions
		ierr = setupPetsc3DMonitor(ts);
		checkPetscError(ierr,
				"PetscSolver::solve: setupPetsc3DMonitor failed.");
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

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Solve the ODE system
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if (ts != NULL && C != NULL) {
		ierr = TSSolve(ts, C);
		checkPetscError(ierr, "PetscSolver::solve: TSSolve failed.");

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		 Write in a file if everything went well or not.
		 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		// Check the option -check_collapse
		PetscBool flagCheck;
		ierr = PetscOptionsHasName(petscOptions, NULL, "-check_collapse",
				&flagCheck);
		checkPetscError(ierr,
				"PetscSolver::solve: PetscOptionsHasName (-check_collapse) failed.");
		if (flagCheck) {
			// Open the output file
			std::ofstream outputFile;
			outputFile.open("solverStatus.txt");

			// Get the converged reason from PETSc
			TSConvergedReason reason;
			ierr = TSGetConvergedReason(ts, &reason);
			checkPetscError(ierr,
					"PetscSolver::solve: TSGetConvergedReason failed.");

			// Write it
			if (reason == TS_CONVERGED_EVENT)
				outputFile << "collapsed" << std::endl;
			else if (reason == TS_DIVERGED_NONLINEAR_SOLVE
					|| reason == TS_DIVERGED_STEP_REJECTED)
				outputFile << "diverged" << std::endl;
			else
				outputFile << "good" << std::endl;

			outputFile.close();
		}
	} else {
		throw std::string(
				"PetscSolver Exception: Unable to solve! Data not configured properly.");
	}

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

	return;
}

void PetscSolver::finalize(bool isStandalone) {
	PetscErrorCode ierr;

	// Switch on the number of dimensions to set the monitors
	int dim = getSolverHandler().getDimension();
	switch (dim) {
	case 0:
		// One dimension
		ierr = reset0DMonitor();
		checkPetscError(ierr, "PetscSolver::finalize: reset0DMonitor failed.");
		break;
	case 1:
		// One dimension
		ierr = reset1DMonitor();
		checkPetscError(ierr, "PetscSolver::finalize: reset1DMonitor failed.");
		break;
	case 2:
		// Two dimensions
		ierr = reset2DMonitor();
		checkPetscError(ierr, "PetscSolver::finalize: reset2DMonitor failed.");
		break;
	case 3:
		// Three dimensions
		ierr = reset3DMonitor();
		checkPetscError(ierr, "PetscSolver::finalize: reset3DMonitor failed.");
		break;
	default:
		throw std::string("PetscSolver Exception: Wrong number of dimensions "
				"to reset the monitors.");
	}

	if (isStandalone) {
		ierr = PetscFinalize();
		checkPetscError(ierr, "PetscSolver::finalize: PetscFinalize failed.");
	}

	return;
}

} /* end namespace xolotlSolver */
