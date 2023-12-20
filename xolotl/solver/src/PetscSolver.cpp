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
#include <xolotl/solver/monitor/PetscMonitor0D.h>
#include <xolotl/solver/monitor/PetscMonitor1D.h>
#include <xolotl/solver/monitor/PetscMonitor2D.h>
#include <xolotl/solver/monitor/PetscMonitor3D.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>

PETSC_EXTERN PetscErrorCode
PetscVFPrintfDefault(FILE*, const char[], va_list);

/*
 C_t =  -D*C_xx + F(C) + R(C)  from Brian Wirth's SciDAC project.

 D*C_xx  - diffusion
 F(C)    - forcing function; impurities being created.
 R(C)    - reaction terms   (production and dissociation)
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

PetscErrorCode
overridePetscVFPrintf(FILE* fd, const char format[], va_list Argp)
{
	PetscFunctionBeginUser;

	PetscErrorCode ierr;
	if (fd == stdout || fd == stderr) {
		constexpr std::size_t BIG = 1024;
		char buff[BIG];
		std::size_t length;
		ierr = PetscVSNPrintf(buff, BIG, format, &length, Argp);
		CHKERRQ(ierr);
		if (fd == stderr) {
			XOLOTL_LOG_ERR << buff;
		}
		else {
			XOLOTL_LOG << buff;
		}
	}
	else {
		ierr = PetscVFPrintfDefault(fd, format, Argp);
		CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}

// PETSc options object
// PetscOptions* petscOptions = nullptr;

// Help message
static char help[] = "Solves C_t =  -D*C_xx + F(C) + R(C) + from "
					 "Brian Wirth's SciDAC project.\n";

void
PetscSolver::setupInitialConditions(DM da, Vec C, DM oldDA, Vec oldC)
{
	// Initialize the concentrations in the solution vector
	this->solverHandler->initializeConcentration(da, C, oldDA, oldC);
}

/* ------------------------------------------------------------------- */

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
RHSFunction(TS ts, PetscReal ftime, Vec C, Vec F, void* ctx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr =
		static_cast<PetscSolver*>(ctx)->rhsFunction(ts, ftime, C, F);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

/*
 Compute the Jacobian entries based on IFunction() and insert them into the
 matrix
 */
PetscErrorCode
RHSJacobian(TS ts, PetscReal ftime, Vec C, Mat A, Mat J, void* ctx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr =
		static_cast<PetscSolver*>(ctx)->rhsJacobian(ts, ftime, C, A, J);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscSolver::PetscSolver(const options::IOptions& options) :
	Solver(options,
		[&options](core::network::IReactionNetwork& network,
			perf::IPerfHandler& perfHandler)
			-> std::shared_ptr<handler::ISolverHandler> {
			switch (options.getDimensionNumber()) {
			case 0:
				return std::make_shared<handler::PetscSolver0DHandler>(
					network, perfHandler, options);
			case 1:
				return std::make_shared<handler::PetscSolver1DHandler>(
					network, perfHandler, options);
			case 2:
				return std::make_shared<handler::PetscSolver2DHandler>(
					network, perfHandler, options);
			case 3:
				return std::make_shared<handler::PetscSolver3DHandler>(
					network, perfHandler, options);
			default:
				// The asked dimension is not good (e.g. -1, 4)
				throw std::runtime_error(
					"\nxolotlFactory: Bad dimension for the solver handler.");
			}
		})
{
	this->setCommandLineOptions(options.getPetscArg());

	rhsFunctionTimer = perfHandler->getTimer("rhsFunctionTimer");
	rhsJacobianTimer = perfHandler->getTimer("rhsJacobianTimer");
	solveTimer = perfHandler->getTimer("solveTimer");
}

PetscSolver::PetscSolver(
	const std::shared_ptr<handler::ISolverHandler>& _solverHandler) :
	Solver(_solverHandler)
{
	rhsFunctionTimer = perfHandler->getTimer("rhsFunctionTimer");
	rhsJacobianTimer = perfHandler->getTimer("rhsJacobianTimer");
	solveTimer = perfHandler->getTimer("solveTimer");
}

PetscSolver::~PetscSolver()
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
PetscSolver::initialize(int loop, double time, DM oldDA, Vec oldC)
{
	PetscErrorCode ierr;

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Initialize program
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if (!isPetscInitialized()) {
		PetscVFPrintf = overridePetscVFPrintf;
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

	// Create the solver context
	this->solverHandler->createSolverContext(da);

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
	ierr = TSSetRHSFunction(ts, nullptr, RHSFunction, this);
	checkPetscError(ierr, "PetscSolver::initialize: TSSetRHSFunction failed.");
	ierr = TSSetRHSJacobian(ts, nullptr, nullptr, RHSJacobian, this);
	checkPetscError(ierr, "PetscSolver::initialize: TSSetRHSJacobian failed.");
	ierr = TSSetSolution(ts, C);
	checkPetscError(ierr, "PetscSolver::initialize: TSSetSolution failed.");

	// Read the times if the information is in the HDF5 file
	auto fileName = this->solverHandler->getNetworkName();
	double deltaTime = 1.0e-12;
	if (!fileName.empty() and loop == 0) {
		io::XFile xfile(fileName);
		auto concGroup = xfile.getGroup<io::XFile::ConcentrationGroup>();
		if (concGroup and concGroup->hasTimesteps()) {
			auto tsGroup = concGroup->getLastTimestepGroup();
			assert(tsGroup);
			std::tie(time, deltaTime) = tsGroup->readTimes();
		}
	}

	// Give the values to the solver
	setCurrentTimes(time, deltaTime);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set solver options
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = TSSetFromOptions(ts);
	checkPetscError(ierr, "PetscSolver::initialize: TSSetFromOptions failed.");

	// Switch on the number of dimensions to set the monitors
	auto dim = this->solverHandler->getDimension();
	switch (dim) {
	case 0:
		// One dimension
		this->monitor =
			std::make_shared<monitor::PetscMonitor0D>(ts, this->solverHandler);
		break;
	case 1:
		// One dimension
		this->monitor =
			std::make_shared<monitor::PetscMonitor1D>(ts, this->solverHandler);
		break;
	case 2:
		// Two dimensions
		this->monitor =
			std::make_shared<monitor::PetscMonitor2D>(ts, this->solverHandler);
		break;
	case 3:
		// Three dimensions
		this->monitor =
			std::make_shared<monitor::PetscMonitor3D>(ts, this->solverHandler);
		break;
	default:
		throw std::runtime_error(
			"PetscSolver Exception: Wrong number of dimensions "
			"to set the monitors.");
	}
	this->monitor->setup(loop);

	// Set the saved data
	if (loop > 0)
		this->monitor->setFlux(
			_nSurf, _nBulk, _previousSurfFlux, _previousBulkFlux);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set initial conditions
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	setupInitialConditions(da, C, oldDA, oldC);

	// Set the output precision for std::out
	std::cout.precision(16);

	// Pop the options
	ierr = PetscOptionsPop();
	checkPetscError(ierr, "PetscSolver::initialize: PetscOptionsPop failed.");

	initTimer->stop();
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
}

std::vector<std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
PetscSolver::getConcVector()
{
	return this->solverHandler->getConcVector(da, C);
}

void
PetscSolver::setConcVector(std::vector<std::vector<
		std::vector<std::vector<std::pair<IdType, double>>>>>& concVector)
{
	return this->solverHandler->setConcVector(da, C, concVector);
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
}

void
PetscSolver::solve()
{
	PetscErrorCode ierr;

	// Initialiaze the converged reason
	TSConvergedReason reason = TS_CONVERGED_USER;
	Vec oldC;
	DM oldDA;
	int loopNumber = 0;
	double time = 0.0;

	// Push the options for the solve
	ierr = PetscOptionsPush(petscOptions);
	checkPetscError(ierr, "PetscSolver::solve: PetscOptionsPush failed.");

	while (reason == TS_CONVERGED_USER) {
		// The interface already initialized the first loop
		if (loopNumber > 0) {
			initialize(loopNumber, time, oldDA, oldC);
		}

		/*
		Solve the ODE system
						- - - - - - - - - - - - - - - - - - - - - - - - - - - -
		- - - - - - */
		if (ts != NULL && C != NULL) {
			// Start the solve Timer
			solveTimer->start();
			if (!petscInitializedHere) {
				// Reset the time step number
				ierr = TSSetStepNumber(ts, 0);
				checkPetscError(
					ierr, "PetscSolver::solve: Reset Step Number failed.");
				// Reset the GB location
				this->solverHandler->initGBLocation(da, C);
			}
			// Start the PETSc Solve
			ierr = TSSolve(ts, C);
			checkPetscError(ierr, "PetscSolver::solve: TSSolve failed.");
			// Stop the timer
			solveTimer->stop();

			// Save some data from the monitors for next loop
			this->monitor->keepFlux(
				_nSurf, _nBulk, _previousSurfFlux, _previousBulkFlux);

			// We are done with the loop
			loopNumber++;

			// Catch the change in surface
			// Get the converged reason from PETSc
			ierr = TSGetConvergedReason(ts, &reason);
			checkPetscError(
				ierr, "PetscSolver::solve: TSGetConvergedReason failed.");
			if (reason == TS_CONVERGED_USER)
				std::cout << "Caught the change of surface!" << std::endl;

			// Save the time
			ierr = TSGetTime(ts, &time);
			checkPetscError(ierr, "PetscSolver::solve: TSGetTime failed.");

			// Save the old DA and associated vector
			PetscInt dof;
			ierr = DMDAGetDof(da, &dof);
			checkPetscError(ierr, "PetscSolver::solve: DMDAGetDof failed.");

			ierr = DMDACreateCompatibleDMDA(da, dof, &oldDA);
			checkPetscError(
				ierr, "PetscSolver::solve: DMDACreateCompatibleDMDA failed.");

			// Save the old vector as a natural one to make the transfer easier
			ierr = DMDACreateNaturalVector(oldDA, &oldC);
			checkPetscError(
				ierr, "PetscSolver::solve: DMDACreateNaturalVector failed.");
			ierr = DMDAGlobalToNaturalBegin(oldDA, C, INSERT_VALUES, oldC);
			ierr = DMDAGlobalToNaturalEnd(oldDA, C, INSERT_VALUES, oldC);
			checkPetscError(
				ierr, "PetscSolver::solve: DMDAGlobalToNatural failed.");
		}
		else {
			throw std::string("PetscSolver Exception: Unable to solve! Data "
							  "not configured properly.");
		}
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			 Write in a file if everything went well or not.
			 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   -*/

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

PetscErrorCode
PetscSolver::rhsFunction(TS ts, PetscReal ftime, Vec C, Vec F)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// Start the RHSFunction Timer
	rhsFunctionTimer->start();

	// Get the local data vector from PETSc
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
	this->solverHandler->updateConcentration(ts, localC, F, ftime);

	// Stop the RHSFunction Timer
	rhsFunctionTimer->stop();

	// Return the local vector
	ierr = DMRestoreLocalVector(da, &localC);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscSolver::rhsJacobian(TS ts, PetscReal ftime, Vec C, Mat A, Mat J)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// Start the RHSJacobian timer
	rhsJacobianTimer->start();

	// Get the matrix from PETSc
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
	this->solverHandler->computeJacobian(ts, localC, J, ftime);

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
	rhsJacobianTimer->stop();

	PetscFunctionReturn(0);
}
} /* end namespace solver */
} /* end namespace xolotl */
