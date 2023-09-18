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

	if (fd == stdout || fd == stderr) {
		constexpr std::size_t BIG = 1024;
		char buff[BIG];
		std::size_t length;
		PetscCall(PetscVSNPrintf(buff, BIG, format, &length, Argp));
		if (fd == stderr) {
			XOLOTL_LOG_ERR << buff;
		}
		else {
			XOLOTL_LOG << buff;
		}
	}
	else {
		PetscCall(PetscVFPrintfDefault(fd, format, Argp));
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
	PetscCall(static_cast<PetscSolver*>(ctx)->rhsFunction(ts, ftime, C, F));
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
	PetscCall(static_cast<PetscSolver*>(ctx)->rhsJacobian(ts, ftime, C, A, J));
	PetscFunctionReturn(0);
}

PetscSolver::PetscSolver(const options::IOptions& options) :
	Solver(options,
		[&options](core::network::IReactionNetwork& network)
			-> std::shared_ptr<handler::ISolverHandler> {
			switch (options.getDimensionNumber()) {
			case 0:
				return std::make_shared<handler::PetscSolver0DHandler>(
					network, options);
			case 1:
				return std::make_shared<handler::PetscSolver1DHandler>(
					network, options);
			case 2:
				return std::make_shared<handler::PetscSolver2DHandler>(
					network, options);
			case 3:
				return std::make_shared<handler::PetscSolver3DHandler>(
					network, options);
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
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Initialize program
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if (!isPetscInitialized()) {
		PetscVFPrintf = overridePetscVFPrintf;
		PetscCallVoid(PetscInitialize(NULL, NULL, NULL, help));
		petscInitializedHere = true;
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create the solver options
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscCallVoid(PetscOptionsCreate(&petscOptions));
	PetscCallVoid(
		PetscOptionsInsertString(petscOptions, optionsString.c_str()));
	PetscCallVoid(PetscOptionsPush(petscOptions));

	// Check the option -snes_mf_operator
	PetscBool flagReduced;
	PetscCallVoid(
		PetscOptionsHasName(NULL, NULL, "-snes_mf_operator", &flagReduced));

	// Create the solver context
	this->solverHandler->createSolverContext(da);
	PetscCallVoid(DMSetVecType(da, VECKOKKOS));
	PetscCallVoid(DMSetMatType(da, MATAIJKOKKOS));

	/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Extract global vector from DMDA to hold solution
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscCallVoid(DMCreateGlobalVector(da, &C));

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create timestepping solver context
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	// Get the MPI communicator
	auto xolotlComm = util::getMPIComm();
	PetscCallVoid(TSCreate(xolotlComm, &ts));
	PetscCallVoid(TSSetType(ts, TSARKIMEX));
	PetscCallVoid(TSARKIMEXSetFullyImplicit(ts, PETSC_TRUE));
	PetscCallVoid(TSSetDM(ts, da));
	PetscCallVoid(TSSetProblemType(ts, TS_NONLINEAR));
	PetscCallVoid(TSSetRHSFunction(ts, nullptr, RHSFunction, this));
	PetscCallVoid(TSSetRHSJacobian(ts, nullptr, nullptr, RHSJacobian, this));
	PetscCallVoid(TSSetSolution(ts, C));

	this->solverHandler->initializeSolverContext(da, ts);

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
	PetscCallVoid(TSSetFromOptions(ts));

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
	PetscCallVoid(PetscOptionsPop());
}

void
PetscSolver::setTimes(double finalTime, double dt)
{
	// Get the default values for the dt
	TSAdapt adapt;
	PetscCallVoid(TSGetAdapt(ts, &adapt));
	PetscReal hmin, hmax;
	PetscCallVoid(TSAdaptGetStepLimits(adapt, &hmin, &hmax));
	// Set the new max value
	PetscCallVoid(TSAdaptSetStepLimits(adapt, hmin, dt));

	// Give the final time value to the solver
	PetscCallVoid(TSSetMaxTime(ts, finalTime));
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
	// Get the value from the solver
	double currentTimeStep;
	PetscCallContinue(TSGetTimeStep(ts, &currentTimeStep));

	return currentTimeStep;
}

void
PetscSolver::setCurrentTimes(double time, double dt)
{
	// Give the values to the solver
	PetscCallVoid(TSSetTime(ts, time));
	PetscCallVoid(TSSetTimeStep(ts, dt));
}

void
PetscSolver::solve()
{
	// Initialiaze the converged reason
	TSConvergedReason reason = TS_CONVERGED_USER;
	Vec oldC;
	DM oldDA;
	int loopNumber = 0;
	double time = 0.0;

	// Push the options for the solve
	PetscCallVoid(PetscOptionsPush(petscOptions));

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
				PetscCallVoid(TSSetStepNumber(ts, 0));
				// Reset the GB location
				this->solverHandler->initGBLocation(da, C);
			}
			// Start the PETSc Solve
			PetscCallVoid(TSSolve(ts, C));
			// Stop the timer
			solveTimer->stop();

			// Save some data from the monitors for next loop
			this->monitor->keepFlux(
				_nSurf, _nBulk, _previousSurfFlux, _previousBulkFlux);

			// We are done with the loop
			loopNumber++;

			// Catch the change in surface
			// Get the converged reason from PETSc
			PetscCallVoid(TSGetConvergedReason(ts, &reason));
			if (reason == TS_CONVERGED_USER)
				std::cout << "Caught the change of surface!" << std::endl;

			// Save the time
			PetscCallVoid(TSGetTime(ts, &time));

			// Save the old DA and associated vector
			PetscInt dof;
			PetscCallVoid(DMDAGetDof(da, &dof));

			PetscCallVoid(DMDACreateCompatibleDMDA(da, dof, &oldDA));

			// Save the old vector as a natural one to make the transfer easier
			PetscCallVoid(DMDACreateNaturalVector(oldDA, &oldC));
			PetscCallVoid(
				DMDAGlobalToNaturalBegin(oldDA, C, INSERT_VALUES, oldC));
			PetscCallVoid(
				DMDAGlobalToNaturalEnd(oldDA, C, INSERT_VALUES, oldC));
		}
		else {
			throw std::string("PetscSolver Exception: Unable to solve! Data "
							  "not configured properly.");
		}
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			 Write in a file if everything went well or not.
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	// Check the option -check_collapse
	PetscBool flagCheck;
	PetscCallVoid(
		PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCheck));
	if (flagCheck) {
		// Open the output file
		std::ofstream outputFile;
		outputFile.open("solverStatus.txt");

		// Get the converged reason from PETSc
		TSConvergedReason reason;
		PetscCallVoid(TSGetConvergedReason(ts, &reason));

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
	PetscCallVoid(PetscOptionsPop());
}

bool
PetscSolver::getConvergenceStatus()
{
	// Get the converged reason from PETSc
	TSConvergedReason reason;
	PetscCallContinue(TSGetConvergedReason(ts, &reason));

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
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Free work space.
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscCallVoid(PetscOptionsDestroy(&petscOptions));
	PetscCallVoid(VecDestroy(&C));
	PetscCallVoid(TSDestroy(&ts));
	PetscCallVoid(DMDestroy(&da));

	if (petscInitializedHere) {
		PetscCallVoid(PetscFinalize());
	}
}

double
PetscSolver::getXolotlTime()
{
	// The most recent time that Xolotl converged
	PetscReal CurrentXolotlTime;
	PetscCallContinue(TSGetTime(ts, &CurrentXolotlTime));
	return CurrentXolotlTime;
}

PetscErrorCode
PetscSolver::rhsFunction(TS ts, PetscReal ftime, Vec C, Vec F)
{
	PetscFunctionBeginUser;

	// Start the RHSFunction Timer
	rhsFunctionTimer->start();

	// Get the local data vector from PETSc
	DM da;
	PetscCall(TSGetDM(ts, &da));
	Vec localC;
	PetscCall(DMGetLocalVector(da, &localC));

	// Scatter ghost points to local vector, using the 2-step process
	// DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
	// By placing code between these two statements, computations can be
	// done while messages are in transition.
	PetscCall(DMGlobalToLocalBegin(da, C, INSERT_VALUES, localC));
	PetscCall(DMGlobalToLocalEnd(da, C, INSERT_VALUES, localC));

	// Set the initial values of F
	PetscCall(VecSet(F, 0.0));

	// Compute the new concentrations
	this->solverHandler->updateConcentration(ts, localC, F, ftime);

	// Stop the RHSFunction Timer
	rhsFunctionTimer->stop();

	// Return the local vector
	PetscCall(DMRestoreLocalVector(da, &localC));

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscSolver::rhsJacobian(TS ts, PetscReal ftime, Vec C, Mat A, Mat J)
{
	PetscFunctionBeginUser;

	// Start the RHSJacobian timer
	rhsJacobianTimer->start();

	// Get the matrix from PETSc
	PetscCall(MatZeroEntries(J));
	DM da;
	PetscCall(TSGetDM(ts, &da));
	Vec localC;
	PetscCall(DMGetLocalVector(da, &localC));

	// Get the complete data array
	PetscCall(DMGlobalToLocalBegin(da, C, INSERT_VALUES, localC));
	PetscCall(DMGlobalToLocalEnd(da, C, INSERT_VALUES, localC));

	// Get the solver handler
	this->solverHandler->computeJacobian(ts, localC, J, ftime);

	// Return the local vector
	PetscCall(DMRestoreLocalVector(da, &localC));

	// PetscCall(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
	// PetscCall(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));

	if (A != J) {
		PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
		PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
	}

	// PetscCall(MatView(J, PETSC_VIEWER_STDOUT_WORLD));

	// Stop the RHSJacobian timer
	rhsJacobianTimer->stop();

	PetscFunctionReturn(0);
}
} /* end namespace solver */
} /* end namespace xolotl */
