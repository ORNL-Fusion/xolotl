// Includes
#include "PetscSolver.h"
#include <xolotlPerf.h>
#include <HDF5NetworkLoader.h>
#include "TemperatureHandler.h"
#include <MathUtils.h>
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <string>
#include <unordered_map>
#include <HDF5Utils.h>

using namespace xolotlCore;

/*
 C_t =  -D*C_xx + F(C) + R(C) + D(C) from Brian Wirth's SciDAC project.

 D*C_xx  - diffusion of He[1-5] and V[1] and I[1]
 F(C)  -   forcing function; He being created.
 R(C)  -   reaction terms   (clusters combining)
 D(C)  -   dissociation terms (cluster breaking up)

 Sample Options:
 -ts_monitor_draw_solution               -- plot the solution for each concentration as a function of x each in a separate 1d graph
 -draw_fields_by_name 1-He-2-V,1-He 	 -- only plot the solution for these two concentrations
 -mymonitor                              -- plot the concentrations of He and V as a function of x and cluster size (2d contour plot)
 -da_refine <n=1,2,...>                  -- run on a finer grid
 -da_grid_x <nx>						 -- number of grid points in the x direction
 -ts_max_steps maxsteps                  -- maximum number of time-steps to take
 -ts_final_time time                     -- maximum time to compute to

 Rules for maximum number of He allowed for V in cluster
 */

namespace xolotlSolver {

//Counter for the number of times RHSFunction is called.
std::shared_ptr<xolotlPerf::IEventCounter> RHSFunctionCounter;

//Counter for the number of times RHSJacobian is called.
std::shared_ptr<xolotlPerf::IEventCounter> RHSJacobianCounter;

//Timer for RHSFunction()
std::shared_ptr<xolotlPerf::ITimer> RHSFunctionTimer;
//Timer for how long it takes to compute ODE terms for each grid point
std::shared_ptr<xolotlPerf::ITimer> computeODEtermPerGP;
//Timer for how long it takes to compute new fluxes in RHSFunction()
//std::shared_ptr<xolotlPerf::ITimer> computeNewFluxes;
//Timer for how long it takes to compute the incident flux in RHSFunction()
//std::shared_ptr<xolotlPerf::ITimer> computeIncidentFluxTimer;
//Timers for how long it takes to compute He, V, and I diffusion in RHSFunction()
//std::shared_ptr<xolotlPerf::ITimer> computeDiffusionTimer;

//Timer for RHSJacobian()
std::shared_ptr<xolotlPerf::ITimer> RHSJacobianTimer;
//Timer for computing reaction term partials per grid point in RHSJacobian()
std::shared_ptr<xolotlPerf::ITimer> computeReactionTermPartials;
//Timer for computing Jacobian diffusion terms per grid point in RHSJacobian()
std::shared_ptr<xolotlPerf::ITimer> computeJacobianDiffusionTerms;
std::shared_ptr<xolotlPerf::ITimer> updateJacobianCol;

//Timer for how long it takes to solve the ODE system in the function solve()
std::shared_ptr<xolotlPerf::ITimer> solveODEsystem;

//! Help message
static char help[] =
		"Solves C_t =  -D*C_xx + F(C) + R(C) + D(C) from Brian Wirth's SciDAC project.\n";

// ----- GLOBAL VARIABLES ----- //

// Allocate the static network
std::shared_ptr<PSIClusterReactionNetwork> PetscSolver::network;
// Allocate the static flux handler
std::shared_ptr<IFluxHandler> PetscSolver::fluxHandler;
// Allocate the static temperature handler
std::shared_ptr<ITemperatureHandler> PetscSolver::temperatureHandler;

extern PetscErrorCode RHSFunction(TS, PetscReal, Vec, Vec, void*);
extern PetscErrorCode RHSJacobian(TS, PetscReal, Vec, Mat, Mat);
extern PetscErrorCode setupPetscMonitor(TS);
extern void computeRetention(TS, Vec);

TS ts; /* nonlinear solver */
Vec C; /* solution */
PetscErrorCode ierr;
DM da; /* manages the grid data */
PetscInt He, *ofill, *dfill;

/**
 * A map for storing the dfill configuration and accelerating the formation of
 * the Jacobian. Its keys are reactant/cluster ids and its values are integer
 * vectors of the column ids that are marked as connected for that cluster in
 * the dfill array.
 */
static std::unordered_map<int, std::vector<int> > dFillMap;

/**
 * A pointer to all of the reactants in the network. It is retrieved from the
 * network after it is set.
 */
static std::shared_ptr<std::vector<Reactant *>> allReactants;

/**
 * The last temperature on the grid. In the future this will have to be an
 * array or map, but for now the temperature is isotropic.
 */
static double lastTemperature = 0.0;

/**
 * A vector for holding the partial derivatives of one cluster. It is sized in
 * the solve() operation.
 *
 * The vector is used for every cluster and immediately reset to zero before
 * being used for the next. This allows the acquisition of the partial
 * derivatives to take up minimal memory and require no additional dynamic
 * allocations.
 */
static std::vector<double> clusterPartials;

/**
 * A vector for holding the partial derivatives for one cluster in the order
 * that PETSc expects. It is sized in the solve() operation.
 *
 * The vector is used for every cluster and immediately reset to zero before
 * being used for the next. This allows the acquisition of the partial
 * derivatives to take up minimal memory and require no additional dynamic
 * allocations.
 */
static std::vector<double> reactingPartialsForCluster;

/* ----- Error Handling Code ----- */

/**
 * This operation checks a Petsc error code and converts it to a bool.
 * @param errorCode The Petsc error code.
 * @return True if everything is OK, false otherwise.
 */
static inline bool checkPetscError(PetscErrorCode errorCode) {
	CHKERRQ(errorCode);
}

/**
 * This operation "returns" in a way that Petsc expects.
 * @return The return code from Petsc.
 */
static inline int petscReturn() {
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "setupInitialConditions"
PetscErrorCode PetscSolver::setupInitialConditions(DM da, Vec C) {

// Local Declarations
	PetscErrorCode ierr;
	PetscInt i, nI, nHe, nV, xs, xm, Mx, cnt = 0;
	PetscScalar *concentrations;
	char string[16];
	int size = allReactants->size();
	double * concOffset;
	std::map<std::string, int> composition;

	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Get the total number of grid points specified by the command line option
	PetscInt numOfxGridPoints;
	PetscBool flg;
	PetscOptionsGetInt(NULL, "-da_grid_x", &numOfxGridPoints, &flg);
	if (!flg)
		numOfxGridPoints = 8.0;

	// Setup some step size variables
	PetscReal hx;
	hx = numOfxGridPoints / (PetscReal) (Mx - 1);
//	// Display the number of grid points that will be used
//	std::cout << "\nNumber of x grid points = " << numOfxGridPoints << std::endl;
//	std::cout << "Number of grid points = " << Mx << std::endl;
//	std::cout << "Step size hx = " << hx << std::endl;
	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = PetscSolver::getFluxHandler();
	fluxHandler->initializeFluxHandler(Mx, hx);

	/* Name each of the concentrations */
	for (i = 0; i < size; i++) {
		composition = allReactants->at(i)->getComposition();
		nHe = composition["He"];
		nV = composition["V"];
		nI = composition["I"];
		ierr = PetscSNPrintf(string, 16, "He-%d,V-%d,I-%d", nHe, nV, nI);
		checkPetscError(ierr);
		ierr = DMDASetFieldName(da, cnt++, string);
		checkPetscError(ierr);
	}

	/*
	 Get pointer to vector data
	 */
	ierr = DMDAVecGetArray(da, C, &concentrations);
	checkPetscError(ierr);

	/*
	 Get local grid boundaries
	 */
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Get the name of the HDF5 file to read the concentrations from
	std::shared_ptr<HDF5NetworkLoader> HDF5Loader
		= std::dynamic_pointer_cast<HDF5NetworkLoader> (networkLoader);
	auto fileName = HDF5Loader->getFilename();

	// Get the last time step written in the HDF5 file
	int tempTimeStep = -2;
	HDF5Utils::hasConcentrationGroup(fileName, tempTimeStep);

	// Loop on all the grid points
	for (i = xs; i < xs + xm; i++) {
		concOffset = concentrations + size * i;
		// Loop on all the clusters to initialize at 0.0
		for (int k = 0; k < size; k++) {
			concOffset[k] = 0.0;
		}

//		// Uncomment this for debugging
//		if (i > 0) {
//			for (int k = 0; k < size; k++) {
//				concOffset[k] = 0.001;
//			}
//		}

		if (tempTimeStep >= 0) {
			// Read the concentrations from the HDF5 file
			auto concVector = HDF5Utils::readGridPoint(fileName,
					tempTimeStep, i);

			// Loop on the concVector size
			for (int k = 0; k < concVector.size(); k++) {
				concOffset[(int) concVector.at(k).at(0)] = concVector.at(k).at(1);
			}
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArray(da, C, &concentrations);
	checkPetscError(ierr);
	PetscFunctionReturn(0);
}

void getIncomingHeFlux(PSICluster * cluster, std::vector<double> gridPos,
		PetscReal curTime, PetscScalar *updatedConcOffset) {

	int reactantIndex = 0;
	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = PetscSolver::getFluxHandler();

	// Get the composition of the cluster
	auto thisComp = cluster->getComposition();
	// Create the composition vector for the cluster
	std::vector<int> compVec = { thisComp["He"], thisComp["V"], thisComp["I"] };
	// Only update the concentration if the cluster exists
	if (cluster) {
		reactantIndex = cluster->getId() - 1;
		// Calculate the incident flux
		auto incidentFlux = fluxHandler->getIncidentFlux(compVec, gridPos,
				curTime);
		// Update the concentration of the cluster
		updatedConcOffset[reactantIndex] += 1.0E4 * PetscMax(0.0, incidentFlux);
		// where incidentFlux = 0.0006 * x * x * x - 0.0087 * x * x + 0.0300 * x
	}

	return;
}

void computeDiffusion(PSICluster * cluster, double temp, PetscReal sx,
		PetscScalar *concOffset, PetscScalar *leftConcOffset,
		PetscScalar *rightConcOffset, PetscScalar *updatedConcOffset) {

	int reactantIndex = 0;
	// Dummy variables to keep the code clean
	double oldConc = 0.0, oldLeftConc = 0.0, oldRightConc = 0.0, conc = 0.0,
			flux = 0.0;

	reactantIndex = cluster->getId() - 1;
	// Get the concentrations
	oldConc = concOffset[reactantIndex];
	oldLeftConc = leftConcOffset[reactantIndex];
	oldRightConc = rightConcOffset[reactantIndex];

	// Use a simple midpoint stencil to compute the concentration
	conc = cluster->getDiffusionCoefficient()
			* (-2.0 * oldConc + oldLeftConc + oldRightConc) * sx;

	// Update the concentration of the cluster
	updatedConcOffset[reactantIndex] += conc;

}

/* ------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
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
PetscErrorCode RHSFunction(TS ts, PetscReal ftime, Vec C, Vec F, void *ptr) {

	// increment the event counter monitoring this function
	RHSFunctionCounter->increment();

	// Important petsc stuff (related to the grid mostly)
	DM da;
	PetscErrorCode ierr;
	PetscInt xi, Mx, xs, xm;
	PetscReal hx, sx, x;
	// Pointers to the Petsc arrays that start at the beginning (xs) of the
	// local array!
	PetscReal *concs, *updatedConcs;
	Vec localC;
	// Loop variables
	int size = 0, reactantIndex = 0;
	// Handy pointers to keep the code clean
	PSICluster * heCluster = NULL, *vCluster = NULL, *iCluster = NULL,
			*cluster = NULL;
	// The following pointers are set to the first position in the conc or
	// updatedConc arrays that correspond to the beginning of the data for the
	// current gridpoint. They are accessed just like regular arrays.
	PetscScalar * concOffset, *leftConcOffset, *rightConcOffset,
			*updatedConcOffset;
	// Dummy variables to keep the code clean
	double flux = 0.0;
	// Get the handle to the network that will be used to compute fluxes.
	auto network = PetscSolver::getNetwork();
	// Some required properties
	auto props = network->getProperties();
	int numHeClusters = std::stoi(props["numHeClusters"]);

	// Get the local data vector from petsc
	PetscFunctionBeginUser;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);
	ierr = DMGetLocalVector(da, &localC);
	checkPetscError(ierr);
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Get the total number of grid points specified by the command line option
	PetscInt numOfxGridPoints;
	PetscBool flg;
	PetscOptionsGetInt(NULL, "-da_grid_x", &numOfxGridPoints, &flg);
	if (!flg)
		numOfxGridPoints = 8.0;

	// Setup some step size variables
	hx = numOfxGridPoints / (PetscReal) (Mx - 1);
	sx = 1.0 / (hx * hx);

	// Scatter ghost points to local vector, using the 2-step process
	// DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
	// By placing code between these two statements, computations can be
	// done while messages are in transition.
	ierr = DMGlobalToLocalBegin(da, C, INSERT_VALUES, localC);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalEnd(da, C, INSERT_VALUES, localC);
	checkPetscError(ierr);

	// Set the initial values of F
	ierr = VecSet(F, 0.0);
	checkPetscError(ierr);

	// Get pointers to vector data
	ierr = DMDAVecGetArray(da, localC, &concs);
	checkPetscError(ierr);
	ierr = DMDAVecGetArray(da, F, &updatedConcs);
	checkPetscError(ierr);

	//Get local grid boundaries
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Variable to represent the real, or current, time
	PetscReal realTime;
	// Get the current time
	ierr = TSGetTime(ts, &realTime);

	// Loop over grid points computing ODE terms for each grid point
	size = network->size();
	computeODEtermPerGP->start();
	for (xi = xs; xi < xs + xm; xi++) {
		x = xi * hx;

		// Vector representing the position at which the flux will be calculated
		// Currently we are only in 1D
		std::vector<double> gridPosition = { x, 0, 0 };

		// Get the temperature handler that will be used to compute fluxes.
		auto temperatureHandler = PetscSolver::getTemperatureHandler();
		auto temperature = temperatureHandler->getTemperature(gridPosition,
				realTime);

		// Update the network if the temperature changed
		if (!xolotlCore::equal(temperature, lastTemperature)) {
			network->setTemperature(temperature);
			lastTemperature = temperature;
		}
				
//		xi = 1; // Uncomment this line for debugging in a single cell.

		// Compute the middle, left, right and new array offsets
		concOffset = concs + size * xi;
		leftConcOffset = concs + size * (xi - 1);
		rightConcOffset = concs + size * (xi + 1);
		updatedConcOffset = updatedConcs + size * xi;

		// Copy data into the PSIClusterReactionNetwork so that it can
		// compute the fluxes properly. The network is only used to compute the
		// fluxes and hold the state data from the last time step. I'm reusing
		// it because it cuts down on memory significantly (about 400MB per
		// grid point) at the expense of being a little tricky to comprehend.
		network->updateConcentrationsFromArray(concOffset);

		// ----- Account for flux of incoming He by computing forcing that
		// produces He of cluster size 1 -----
		// Crude cubic approximation of graph from Tibo's notes
		heCluster = (PSICluster *) network->get("He", 1);
		//computeIncidentFluxTimer->start();
		getIncomingHeFlux(heCluster, gridPosition, realTime, updatedConcOffset);
		//computeIncidentFluxTimer->stop();

		// ---- Compute diffusion over the locally owned part of the grid -----
		//computeDiffusionTimer->start();

		// He clusters larger than 5 do not diffuse -- they are immobile
		for (int i = 1; i < PetscMin(numHeClusters + 1, 6); i++) {
			// Get the reactant index
			heCluster = (PSICluster *) network->get("He", i);

			// Only update the concentration if the cluster exists
			if (heCluster) {
				computeDiffusion(heCluster, temperature, sx, concOffset,
						leftConcOffset, rightConcOffset, updatedConcOffset);
			}
		}

		// ----- Vacancy Diffusion -----
		// Only vacancy clusters of size 1 diffuse, so grab 1V.
		vCluster = (PSICluster *) network->get("V", 1);
		// Only update the concentration if the cluster exists
		if (vCluster) {
			computeDiffusion(vCluster, temperature, sx, concOffset,
					leftConcOffset, rightConcOffset, updatedConcOffset);
		}

		// ----- Interstitial Diffusion -----
		// Get 1I from the new network and gets its position in the array
		iCluster = (PSICluster *) network->get("I", 1);
		// Only update the concentration if the clusters exist
		if (iCluster) {
			computeDiffusion(iCluster, temperature, sx, concOffset,
					leftConcOffset, rightConcOffset, updatedConcOffset);
		}
		//computeDiffusionTimer->stop();

		// ----- Compute all of the new fluxes -----
		//computeNewFluxes->start();
		for (int i = 0; i < size; i++) {
			cluster = (PSICluster *) allReactants->at(i);
			// Compute the flux
			flux = cluster->getTotalFlux(temperature);
			// Update the concentration of the cluster
			reactantIndex = cluster->getId() - 1;
			updatedConcOffset[reactantIndex] += flux;
//			std::cout << "New flux = " << flux << " "
//					<< cluster->getConcentration() << std::endl;
		}
		//computeNewFluxes->stop();

		// Boundary conditions
		if (xi == 0) {
			for (int i = 0; i < size; i++) {
				updatedConcOffset[i] = 1.0 * concs[i];
			}
		}

//		for (int i = 0; i < size; i++) {
//			std::cout << updatedConcOffset[i] << std::endl;
//		}

		// Uncomment this line for debugging in a single cell.
//		break;
	}
	computeODEtermPerGP->stop();

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArray(da, localC, &concs);
	checkPetscError(ierr);
	ierr = DMDAVecRestoreArray(da, F, &updatedConcs);
	checkPetscError(ierr);
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr);
	PetscFunctionReturn(0);

}

PetscErrorCode callRHSFunction(TS ts, PetscReal ftime, Vec C, Vec F,
		void *ptr) {
	PetscErrorCode ierr;
	RHSFunctionTimer->start();
	ierr = RHSFunction(ts, ftime, C, F, &ptr);
	RHSFunctionTimer->stop();

	return ierr;
}

void computePartialsForDiffusion(PSICluster * cluster, double temp,
		PetscReal sx, PetscReal val[6], PetscInt row[3], PetscInt col[3],
		PetscInt xi, PetscInt xs, int size) {

	int reactantIndex = 0;
	double diffCoeff = 0.0;

	diffCoeff = cluster->getDiffusionCoefficient();
	// Compute the partial derivatives for diffusion of this cluster
	val[0] = diffCoeff * sx;
	val[1] = -2.0 * diffCoeff * sx;
	val[2] = diffCoeff * sx;
	// Get the reactant index
	reactantIndex = cluster->getId() - 1;
	// Set the row and column indices. These indices are computed
	// by using xi, xi-1 and xi+1 and the arrays are shifted to
	// (xs+1)*size to properly account for the neighboring ghost
	// cells.
	row[0] = (xi - xs + 1) * size + reactantIndex;
	col[0] = ((xi - 1) - xs + 1) * size + reactantIndex;
	col[1] = (xi - xs + 1) * size + reactantIndex;
	col[2] = ((xi + 1 + 1) - xs) * size + reactantIndex;

	// Boundary conditions
	if (xi == 0) {
		val[0] = 0.0;
		val[1] = 0.0;
		val[2] = 0.0;
	}
}

#undef __FUNCT__
#define __FUNCT__ "RHSJacobian"
/*
 Compute the Jacobian entries based on IFuction() and insert them into the matrix
 */
PetscErrorCode RHSJacobian(TS ts, PetscReal ftime, Vec C, Mat A, Mat J,
		void *ptr) {

	// increment the event counter monitoring this function
	RHSJacobianCounter->increment();

	DM da;
	PetscErrorCode ierr;
	PetscInt xi, Mx, xs, xm, i;
	PetscInt row[3], col[3];
	PetscReal hx, sx, val[6];
	PetscReal *concs, *updatedConcs;
	double * concOffset;
	Vec localC;
	static PetscBool initialized = PETSC_FALSE;
	// Get the network
	auto network = PetscSolver::getNetwork();
	// Get the properties
	auto props = network->getProperties();
	int numHeClusters = std::stoi(props["numHeClusters"]);
	int reactantIndex = 0;
	int size = 0;

	// Get the matrix from PETSc
	PetscFunctionBeginUser;
	ierr = MatZeroEntries(J);
	checkPetscError(ierr);
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);
	ierr = DMGetLocalVector(da, &localC);
	checkPetscError(ierr);
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Get the total number of grid points specified by the command line option
	PetscInt numOfxGridPoints;
	PetscBool flg;
	PetscOptionsGetInt(NULL, "-da_grid_x", &numOfxGridPoints, &flg);
	if (!flg)
		numOfxGridPoints = 8.0;

	// Setup some step size variables
	hx = numOfxGridPoints / (PetscReal) (Mx - 1);
	sx = 1.0 / (hx * hx);

	// Get the complete data array
	ierr = DMGlobalToLocalBegin(da, C, INSERT_VALUES, localC);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalEnd(da, C, INSERT_VALUES, localC);
	checkPetscError(ierr);

	/*
	 The f[] is a dummy, values are never set into it. It is only used to
	 determine the local row for the entries in the Jacobian
	 */
	ierr = DMDAVecGetArray(da, localC, &concs);
	checkPetscError(ierr);
	ierr = DMDAVecGetArray(da, localC, &updatedConcs);
	checkPetscError(ierr);
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Store network size for both the linear and nonlinear parts of the
	// computation.
	size = network->size();

	// Variable to represent the real, or current, time
	PetscReal realTime;
	// Get the current time
	ierr = TSGetTime(ts, &realTime);

	// Only compute the linear part of the Jacobian once
	if (!initialized) {

		computeJacobianDiffusionTerms->start();
		/*
		 Loop over grid points computing Jacobian terms for diffusion at each
		 grid point
		 */
		for (xi = xs; xi < xs + xm; xi++) {
			// Vector representing the position at which the flux will be calculated
			// Currently we are only in 1D
			std::vector<double> gridPosition = { xi * hx, 0, 0 };
			// Get the temperature handler that will be used to compute fluxes.
			auto temperatureHandler = PetscSolver::getTemperatureHandler();
			auto temperature = temperatureHandler->getTemperature(gridPosition,
					realTime);

//			xi = 1; // Uncomment this line for debugging in a single cell

			// Copy data into the PSIClusterReactionNetwork so that it can
			// compute the new concentrations.
			concOffset = concs + size * xi;
			network->updateConcentrationsFromArray(concOffset);

			/* -------------------------------------------------------------
			 ---- Compute diffusion over the locally owned part of the grid
			 */

			/* He clusters larger than 5 do not diffuse -- they are immobile */
			// ---- Compute diffusion over the locally owned part of the grid -----
			for (i = 1; i < PetscMin(numHeClusters + 1, 6); i++) {
				// Get the cluster
				auto psiCluster = (PSICluster *) network->get("He", i);
				computePartialsForDiffusion(psiCluster, temperature, sx, val,
						row, col, xi, xs, size);
				ierr = MatSetValuesLocal(J, 1, row, 3, col, val, ADD_VALUES);
				checkPetscError(ierr);
			}

			/* 1V and 1I are the only other clusters that diffuse */

			// ----- Vacancy Diffusion -----
			auto psiCluster = (PSICluster *) network->get("V", 1);
			if (psiCluster) {
				computePartialsForDiffusion(psiCluster, temperature, sx, val,
						row, col, xi, xs, size);
				ierr = MatSetValuesLocal(J, 1, row, 3, col, val, ADD_VALUES);
				checkPetscError(ierr);
			}

			// ----- Interstitial Diffusion -----
			psiCluster = (PSICluster *) network->get("I", 1);
			if (psiCluster) {
				computePartialsForDiffusion(psiCluster, temperature, sx, val,
						row, col, xi, xs, size);
				ierr = MatSetValuesLocal(J, 1, row, 3, col, val, ADD_VALUES);
				checkPetscError(ierr);
			}

//			break;   // Uncomment this line for debugging in a single cell.
		}
		computeJacobianDiffusionTerms->stop();

		ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
		checkPetscError(ierr);
		ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
		checkPetscError(ierr);
//		ierr = MatSetOption(J, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
//		checkPetscError(ierr);
		ierr = MatStoreValues(J);
		checkPetscError(ierr);
//		MatSetFromOptions(J);
		initialized = PETSC_TRUE;
		// Debug line for viewing the matrix
		//MatView(J, PETSC_VIEWER_STDOUT_WORLD);
	} else {
		ierr = MatRetrieveValues(J);
		checkPetscError(ierr);
	}

	/* ----- Compute the partial derivatives for the reaction term at each
	 * grid point ----- */

	// Create a new row array of size n
	PetscInt localPDColIds[size];
	PetscInt rowId = 0;
	int pdColIdsVectorSize = 0;

	computeReactionTermPartials->start();
	// Loop over the grid points
	for (xi = xs; xi < xs + xm; xi++) {
		// Vector representing the position at which the flux will be calculated
		// Currently we are only in 1D
		std::vector<double> gridPosition = { xi * hx, 0, 0 };
		// Get the temperature handler that will be used to compute fluxes.
		auto temperatureHandler = PetscSolver::getTemperatureHandler();
		auto temperature = temperatureHandler->getTemperature(gridPosition,
				realTime);

//		xi = 1; // Uncomment this line for debugging in a single cell

		// Copy data into the PSIClusterReactionNetwork so that it can
		// compute the new concentrations.
		concOffset = concs + size * xi;
		network->updateConcentrationsFromArray(concOffset);
		updateJacobianCol->start();
		// Update the column in the Jacobian that represents each reactant
		for (int i = 0; i < size; i++) {
			auto reactant = allReactants->at(i);
			// Get the reactant index
			reactantIndex = reactant->getId() - 1;
			// Get the column id
			rowId = (xi - xs + 1) * size + reactantIndex;
			// Get the partial derivatives
			reactant->getPartialDerivatives(temperature, clusterPartials);
			// Get the list of column ids from the map
			auto pdColIdsVector = dFillMap.at(reactantIndex);
			//Number of partial derivatives
			pdColIdsVectorSize = pdColIdsVector.size();
			// Loop over the list of column ids
			for (int j = 0; j < pdColIdsVectorSize; j++) {
				// Calculate the appropriate index to match the dfill array
				// configuration
				localPDColIds[j] = (xi - xs + 1) * size + pdColIdsVector[j];
				// Get the partial derivative from the array of all of the partials
				reactingPartialsForCluster[j] =
						clusterPartials[pdColIdsVector[j]];
				// Reset the cluster partial value to zero. This is much faster
				// than using memset.
				clusterPartials[pdColIdsVector[j]] = 0.0;
			}
			// Update the matrix
			ierr = MatSetValuesLocal(J, 1, &rowId, pdColIdsVectorSize,
					localPDColIds, reactingPartialsForCluster.data(),
					ADD_VALUES);
			checkPetscError(ierr);
		}
		updateJacobianCol->stop();
		// Uncomment this line for debugging in a single cell.
//		break;
	}
	computeReactionTermPartials->stop();

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArray(da, C, &concs);
	checkPetscError(ierr);
	ierr = DMDAVecRestoreArray(da, C, &updatedConcs);
	checkPetscError(ierr);
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr);
	ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
	checkPetscError(ierr);
	ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
	checkPetscError(ierr);

	// Enforce the Boundary conditions. Loop over the grid points and set the
	// conditions.
	for (xi = xs; xi < xs + xm; xi++) {
		if (xi == 0) {
			// Loop on the reactants
			for (int i = 0; i < size; i++) {
				auto reactant = allReactants->at(i);
				// Get the reactant index
				reactantIndex = reactant->getId() - 1;
				// Get the row id
				rowId = (xi - xs + 1) * size + reactantIndex;

				// Get the list of column ids from the map
				auto pdColIdsVector = dFillMap.at(reactantIndex);
				pdColIdsVectorSize = pdColIdsVector.size(); //Number of partial derivatives
				// Loop over the list of column ids
				for (int j = 0; j < pdColIdsVectorSize; j++) {
					// Calculate the appropriate index to match the dfill array configuration
					localPDColIds[j] = (xi - xs + 1) * size + pdColIdsVector[j];
					// Get the partial derivative from the array of all of the partials
					reactingPartialsForCluster[j] = 0.0;
				}

				// Update the matrix
				ierr = MatSetValuesLocal(J, 1, &rowId, pdColIdsVectorSize,
						localPDColIds, reactingPartialsForCluster.data(),
						INSERT_VALUES);
				checkPetscError(ierr);
			}
		}
	}

	// Assemble again
	ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
	checkPetscError(ierr);
	ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
	checkPetscError(ierr);

	if (A != J) {
		ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
		checkPetscError(ierr);
		ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
		checkPetscError(ierr);
	}
	PetscFunctionReturn(0);

}

PetscErrorCode callRHSJacobian(TS ts, PetscReal ftime, Vec C, Mat A, Mat J,
		void *ptr) {
	PetscErrorCode ierr;
	RHSJacobianTimer->start();
	ierr = RHSJacobian(ts, ftime, C, A, J, &ptr);
	RHSJacobianTimer->stop();

	return ierr;
}

#undef __FUNCT__
#define __FUNCT__ "getDiagonalFill"

PetscErrorCode PetscSolver::getDiagonalFill(PetscInt *diagFill,
		int diagFillSize) {

	// Local Declarations
	int i = 0, j = 0, numReactants = network->size(), index = 0, id = 0,
			connectivityLength = 0, size = numReactants * numReactants;
	std::vector<int> connectivity;

	// Fill the diagonal block if the sizes match up
	if (diagFillSize == size) {
		// Get the connectivity for each reactant
		for (i = 0; i < numReactants; i++) {
			// Get the reactant and its connectivity
			auto reactant = allReactants->at(i);
			connectivity = reactant->getConnectivity();
			connectivityLength = connectivity.size();
			// Get the reactant id so that the connectivity can be lined up in
			// the proper column
			id = reactant->getId() - 1;
			// Create the vector that will be inserted into the dFill map
			std::vector<int> columnIds;
			// Add it to the diagonal fill block
			for (j = 0; j < connectivityLength; j++) {
				// The id starts at j*connectivity length and is always offset
				// by the id, which denotes the exact column.
				index = id * connectivityLength + j;
				diagFill[index] = connectivity[j];
				// Add a column id if the connectivity is equal to 1.
				if (connectivity[j] == 1) {
					columnIds.push_back(j);
				}
			}
			// Update the map
			dFillMap[id] = columnIds;
		}
		// Debug output
//		std::cout << "Number of degrees of freedom = " << numReactants
//				<< std::endl;
//		printf("\n");
//		for (i = 0; i < numReactants; i++) {
//			for (j = 0; j < numReactants; j++) {
//				printf("%d ", dfill[i * numReactants + j]);
//			}
//			printf("\n");
//		}
//		printf("\n");
	} else {
		std::string err =
				"PetscSolver Exception: Invalid diagonal block size!\n";
		throw std::string(err);
	}

	return 0;
}

//! The Constructor
PetscSolver::PetscSolver() {
	numCLIArgs = 0;
	CLIArgs = NULL;
}

PetscSolver::PetscSolver(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		handlerRegistry(registry) {

	numCLIArgs = 0;
	CLIArgs = NULL;

	RHSFunctionCounter = handlerRegistry->getEventCounter(
			"Petsc_RHSFunction_Counter");
	RHSJacobianCounter = handlerRegistry->getEventCounter(
			"Petsc_RHSJacobian_Counter");

	RHSFunctionTimer = handlerRegistry->getTimer("RHSFunctionTimer");
	computeODEtermPerGP = handlerRegistry->getTimer("computeODEtermPerGP");
	//computeNewFluxes = handlerRegistry->getTimer("computeNewFluxes");
	//computeIncidentFluxTimer = handlerRegistry->getTimer("computeIncidentFluxTimer");
	//computeDiffusionTimer = handlerRegistry->getTimer("computeDiffusionTimer");
	RHSJacobianTimer = handlerRegistry->getTimer("RHSJacobianTimer");
	computeReactionTermPartials = handlerRegistry->getTimer(
			"computeReactionTermPartials");
	computeJacobianDiffusionTerms = handlerRegistry->getTimer(
			"computeJacobianDiffusionTerms");
	updateJacobianCol = handlerRegistry->getTimer("updateJacobianCol");
	solveODEsystem = handlerRegistry->getTimer("solveODEsystem");

}

//! The Destructor
PetscSolver::~PetscSolver() {

    // std::cerr << "Destroying a PetscSolver" << std::endl;

    // Break "pointer" cycles so that network, clusters, reactants
    // will deallocate when the std::shared_ptrs owning them 
    // are destroyed.
    network->askReactantsToReleaseNetwork();
}

/**
 * This operation transfers the input arguments passed to the program on
 * startup to the solver. These options are static options specified at
 * the start of the program whereas the options passed to setOptions() may
 * change.
 * @param argc The number of command line arguments
 * @param argv The array of command line arguments
 */
void PetscSolver::setCommandLineOptions(int argc, char **argv) {

	numCLIArgs = argc;
	CLIArgs = argv;
}

/**
 * This operation sets the PSIClusterNetworkLoader that should be used by
 * the ISolver to load the ReactionNetwork.
 * @param networkLoader The PSIClusterNetworkLoader that will load the
 * network.
 */
void PetscSolver::setNetworkLoader(
		std::shared_ptr<PSIClusterNetworkLoader> networkLoader) {

	// Store the loader and load the network
	this->networkLoader = networkLoader;
	network = networkLoader->load();

	// Get all of the reactants
	allReactants = network->getAll();

	// Debug output
	// Get the processor id
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);
	if (procId == 1) {
		std::cout << "PETScSolver Message: " << "Master loaded network of size "
				<< network->size() << "." << std::endl;
	}

	return;
}

/**
 * This operation sets the run-time options of the solver. The map is a set
 * of key-value std::string pairs that are interpreted by the solver. These
 * options may change during execution, but it is up to Solvers to monitor
 * the map for changes and they may do so at their discretion.
 * @param options The set of options as key-value pairs with option names
 * for keys and associated values mapped to those keys. A relevant example
 * is "startTime" and "0.01" where both are of type std::string.
 */
void PetscSolver::setOptions(std::map<std::string, std::string> options) {
}

/**
 * This operation sets up the mesh that will be used by the solver and
 * initializes the data on that mesh. This operation will throw an exception
 * of type std::string if the mesh can not be setup.
 */
void PetscSolver::setupMesh() {
}

#undef __FUNCT__
#define __FUNCT__ "initialize"
/**
 * This operation performs all necessary initialization for the solver
 * possibly including but not limited to setting up MPI and loading initial
 * conditions. If the solver can not be initialized, this operation will
 * throw an exception of type std::string.
 */
void PetscSolver::initialize() {
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Initialize program
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscInitialize(&numCLIArgs, &CLIArgs, (char*) 0, help);

	return;
}

#undef __FUNCT__
#define __FUNCT__ "solve"
/**
 * This operation directs the Solver to perform the solve. If the solve
 * fails, it will throw an exception of type std::string.
 */
void PetscSolver::solve(std::shared_ptr<IFluxHandler> fluxHandler,
		std::shared_ptr<ITemperatureHandler> temperatureHandler) {

	// Set the flux handler
	PetscSolver::fluxHandler = fluxHandler;

	// Set the temperature handler
	PetscSolver::temperatureHandler = temperatureHandler;

	// Get the properties
	auto props = network->getProperties();
	int numHeClusters = std::stoi(props["numHeClusters"]);
	// The degrees of freedom should be equal to the number of reactants.
	int dof = network->size();

	// Set the size of the partial derivatives vectors
	clusterPartials.resize(dof, 0.0);
	reactingPartialsForCluster.resize(dof, 0.0);

	// Check the network before getting busy.
	if (!network) {
		throw std::string("PetscSolver Exception: Network not set!");
	}

	// Set the output precision for std::out
	std::cout.precision(16);

	PetscFunctionBeginUser;

	// Get the name of the HDF5 file to read the concentrations from
	std::shared_ptr<HDF5NetworkLoader> HDF5Loader
		= std::dynamic_pointer_cast<HDF5NetworkLoader> (networkLoader);
	auto fileName = HDF5Loader->getFilename();

	// Get starting conditions from HDF5 file
	int gridLength = 0;
	double time = 0.0, deltaTime = 1.0e-8;
	int tempTimeStep = -2;
	HDF5Utils::readHeader(fileName, gridLength);

	// Read the times if the information is in the HDF5 file
	if (HDF5Utils::hasConcentrationGroup(fileName, tempTimeStep)) {
		HDF5Utils::readTimes(fileName, tempTimeStep, time, deltaTime);
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_MIRROR, -8, dof, 1,
	NULL, &da);
	checkPetscError(ierr);

	/* The only spatial coupling in the Jacobian (diffusion) is for the first 5 He, the first V, and the first I.
	 The ofill (thought of as a dof by dof 2d (row-oriented) array represents the nonzero coupling between degrees
	 of freedom at one point with degrees of freedom on the adjacent point to the left or right. A 1 at i,j in the
	 ofill array indicates that the degree of freedom i at a point is coupled to degree of freedom j at the
	 adjacent point. In this case ofill has only a few diagonal entries since the only spatial coupling is regular diffusion. */
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &ofill);
	checkPetscError(ierr);
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &dfill);
	checkPetscError(ierr);
	ierr = PetscMemzero(ofill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr);
	ierr = PetscMemzero(dfill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr);

	// Fill ofill, the matrix of "off-diagonal" elements that represents diffusion, with for He.
	int reactantIndex = 0;
	Reactant * reactant;
	for (int numHe = 1; numHe < PetscMin(numHeClusters + 1, 6); numHe++) {
		reactant = network->get("He", numHe);
		// Only couple if the reactant exists
		if (reactant) {
			// Subtract one from the id to get a unique index between 0 and network->size() - 1
			reactantIndex = reactant->getId() - 1;
			ofill[reactantIndex * dof + reactantIndex] = 1;
		}
	}
	// Now for single V
	reactant = network->get("V", 1);
	// Only couple if the reactant exists
	if (reactant) {
		// Subtract one from the id to get a unique index between 0 and network->size() - 1
		reactantIndex = reactant->getId() - 1;
		ofill[reactantIndex * dof + reactantIndex] = 1;
	}
	// Now for single I
	reactant = network->get("I", 1);
	// Only couple if the reactant exists
	if (reactant) {
		// Subtract one from the id to get a unique index between 0 and network->size() - 1
		reactantIndex = reactant->getId() - 1;
		ofill[reactantIndex * dof + reactantIndex] = 1;
	}

	// Get the diagonal fill
	ierr = getDiagonalFill(dfill, dof * dof);
	checkPetscError(ierr);

	// Load up the block fills
	ierr = DMDASetBlockFills(da, dfill, ofill);

	checkPetscError(ierr);
	// Free the temporary fill arrays
	ierr = PetscFree(ofill);
	checkPetscError(ierr);
	ierr = PetscFree(dfill);
	checkPetscError(ierr);

	/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Extract global vector from DMDA to hold solution
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = DMCreateGlobalVector(da, &C);
	checkPetscError(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create timestepping solver context
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = TSCreate(PETSC_COMM_WORLD, &ts);
	checkPetscError(ierr);
	ierr = TSSetType(ts, TSARKIMEX);
	checkPetscError(ierr);
	ierr = TSARKIMEXSetFullyImplicit(ts, PETSC_TRUE);
	checkPetscError(ierr);
	ierr = TSSetDM(ts, da);
	checkPetscError(ierr);
	ierr = TSSetProblemType(ts, TS_NONLINEAR);
	checkPetscError(ierr);
	ierr = TSSetRHSFunction(ts, NULL, callRHSFunction, NULL);
//ierr = TSSetRHSFunction(ts, NULL, RHSFunction, NULL);
	checkPetscError(ierr);
	ierr = TSSetRHSJacobian(ts, NULL, NULL, callRHSJacobian, NULL);
//ierr = TSSetRHSJacobian(ts, NULL, NULL, RHSJacobian, NULL);
	checkPetscError(ierr);
	ierr = TSSetSolution(ts, C);
	checkPetscError(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set solver options
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = TSSetInitialTimeStep(ts, time, deltaTime);
	checkPetscError(ierr);
//	ierr = TSSetDuration(ts, 100, 50.0);
//	checkPetscError(ierr);
	ierr = TSSetFromOptions(ts);
	checkPetscError(ierr);
	ierr = setupPetscMonitor(ts);
	checkPetscError(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set initial conditions
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = setupInitialConditions(da, C);
	checkPetscError(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Solve the ODE system
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	// Time how long it takes to solve the ODE system
	solveODEsystem->start();  // start the timer
	if (ts != NULL && C != NULL) {
		ierr = TSSolve(ts, C);
		checkPetscError(ierr);

		solveODEsystem->stop();  // stop the timer

		// Flags to launch the monitors or not
		PetscBool flagRetention;

		// Check the option -helium_retention
		ierr = PetscOptionsHasName(NULL, "-helium_retention", &flagRetention);
		checkPetscError(ierr);

		if (flagRetention)
			computeRetention(ts, C);
	} else {
		throw std::string(
				"PetscSolver Exception: Unable to solve! Data not configured properly.");
	}
}

/**
 * This operation performs all necessary finalization for the solver
 * including but not limited to cleaning up memory, finalizing MPI and
 * printing diagnostic information. If the solver can not be finalized,
 * this operation will throw an exception of type std::string.
 */
void PetscSolver::finalize() {

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Free work space.
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = VecDestroy(&C);
	checkPetscError(ierr);
	ierr = TSDestroy(&ts);
	checkPetscError(ierr);
	ierr = DMDestroy(&da);
	checkPetscError(ierr);
	ierr = PetscFinalize();
	if (petscReturn() != 0) {
		throw std::string("PetscSolver Exception: Unable to finalize solve!");
	}

}

} /* end namespace xolotlSolver */
