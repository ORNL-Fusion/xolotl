// Includes
#include "PetscSolver.h"
#include "../xolotlPerf/HandlerRegistryFactory.h"
#include "FitFluxHandler.h"
#include "TemperatureHandler.h"
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <string>
#include <unordered_map>

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
std::shared_ptr<xolotlPerf::ITimer> computeNewFluxes;
//Timer for how long it takes to compute the incident flux in RHSFunction()
std::shared_ptr<xolotlPerf::ITimer> computeIncidentFlux;
//Timers for how long it takes to compute He, V, and I diffusion in RHSFunction()
std::shared_ptr<xolotlPerf::ITimer> computeAllDiffusion;

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
extern PetscErrorCode RHSJacobian(TS, PetscReal, Vec, Mat*, Mat*, MatStructure*,
		void*);
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
	auto reactants = network->getAll();
	int size = reactants->size();
	double * concOffset;
	std::map<std::string, int> composition;

	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr);

	/* Name each of the concentrations */
	for (i = 0; i < size; i++) {
		composition = reactants->at(i)->getComposition();
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

	/*
	 Compute function over the locally owned part of the grid
	 */
	for (i = xs; i < xs + xm; i++) {

		// For boundary conditions all the concentrations are 0
		// at i == 0. Everywhere else, only I1 has a non-zero concentration.
		if (i != 0) {
			// Set the default interstitial concentrations
			auto reactant = network->get("I",1);
			reactant->setConcentration(0.0023);
		}

		// Update the PETSc concentrations array
		concOffset = concentrations + size * i;
		network->fillConcentrationsArray(concOffset);
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArray(da, C, &concentrations);
	checkPetscError(ierr);
	PetscFunctionReturn(0);
}

void incomingHeFlux(std::shared_ptr<PSICluster> cluster,
		std::vector<double> gridPos, PetscReal curTime,
		PetscScalar *updatedConcOffset) {

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
}

void computeDiffusion(std::shared_ptr<PSICluster> cluster, double temp,
		PetscReal sx, PetscScalar *concOffset, PetscScalar *leftConcOffset,
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
	conc = cluster->getDiffusionCoefficient(temp)
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
	std::shared_ptr<PSICluster> heCluster, vCluster, iCluster, cluster;
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
	hx = numOfxGridPoints / (PetscReal)(Mx - 1);
	// Display the number of grid points that will be used
//	std::cout << "\nNumber of x grid points = " << numOfxGridPoints << std::endl;
//	std::cout << "Number of grid points = " << Mx << std::endl;
//	std::cout << "Step size hx = " << hx << std::endl;
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
		std::vector<double> gridPosition = { 0, x, 0 };

		// Get the temperature handler that will be used to compute fluxes.
		auto temperatureHandler = PetscSolver::getTemperatureHandler();
		auto temperature = temperatureHandler->getTemperature(gridPosition,
				realTime);

		//xi = 4; // Debugging

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
		heCluster = std::dynamic_pointer_cast < PSICluster
				> (network->get("He", 1));
		computeIncidentFlux->start();
		incomingHeFlux(heCluster, gridPosition, realTime, updatedConcOffset);
		computeIncidentFlux->stop();

		// ---- Compute diffusion over the locally owned part of the grid -----
		computeAllDiffusion->start();

		// He clusters larger than 5 do not diffuse -- they are immobile
		for (int i = 1; i < PetscMin(numHeClusters + 1, 6); i++) {
			// Get the reactant index
			heCluster = std::dynamic_pointer_cast < PSICluster
					> (network->get("He", i));

			// Only update the concentration if the cluster exists
			if (heCluster) {
				computeDiffusion(heCluster, temperature, sx, concOffset,
						leftConcOffset, rightConcOffset, updatedConcOffset);
			}
		}

		// ----- Vacancy Diffusion -----
		// Only vacancy clusters of size 1 diffuse, so grab 1V.
		vCluster = std::dynamic_pointer_cast < PSICluster
				> (network->get("V", 1));
		// Only update the concentration if the cluster exists
		if (vCluster) {
			computeDiffusion(vCluster, temperature, sx, concOffset,
					leftConcOffset, rightConcOffset, updatedConcOffset);
		}

		// ----- Interstitial Diffusion -----
		// Get 1I from the new network and gets its position in the array
		iCluster = std::dynamic_pointer_cast < PSICluster
				> (network->get("I", 1));
		// Only update the concentration if the clusters exist
		if (iCluster) {
			computeDiffusion(iCluster, temperature, sx, concOffset,
					leftConcOffset, rightConcOffset, updatedConcOffset);
		}
		computeAllDiffusion->stop();

		// ----- Compute all of the new fluxes -----
		computeNewFluxes->start();
		auto reactants = network->getAll();
		for (int i = 0; i < size; i++) {
			cluster = std::dynamic_pointer_cast < PSICluster
					> (reactants->at(i));
			// Compute the flux
			flux = cluster->getTotalFlux(temperature);
			// Update the concentration of the cluster
			reactantIndex = cluster->getId() - 1;
			updatedConcOffset[reactantIndex] += flux;
//			std::cout << "New flux = " << flux << " "
//					<< cluster->getConcentration() << std::endl;
		}
		computeNewFluxes->stop();

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
		//break;
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

void computePartialsForDiffusion(std::shared_ptr<PSICluster> cluster,
		double temp, PetscReal sx, PetscReal val[6], PetscInt row[3],
		PetscInt col[3], PetscInt xi, PetscInt xs, int size) {

	int reactantIndex = 0;
	double diffCoeff = 0.0;

	diffCoeff = cluster->getDiffusionCoefficient(temp);
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
PetscErrorCode RHSJacobian(TS ts, PetscReal ftime, Vec C, Mat *A, Mat *J,
		MatStructure *str, void *ptr) {

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
	std::shared_ptr<PSICluster> psiCluster;
	std::shared_ptr<Reactant> heCluster, reactant;
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>>reactants;
	// Get the network
	auto network = PetscSolver::getNetwork();
	// Get the properties
	auto props = network->getProperties();
	int numHeClusters = std::stoi(props["numHeClusters"]);
	int reactantIndex = 0;
	int size = 0;

	// Get the matrix from PETSc
	PetscFunctionBeginUser;
	ierr = MatZeroEntries(*J);
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
	hx = numOfxGridPoints / (PetscReal)(Mx - 1);
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
			std::vector<double> gridPosition = { 0, xi * hx, 0 };
			// Get the temperature handler that will be used to compute fluxes.
			auto temperatureHandler = PetscSolver::getTemperatureHandler();
			auto temperature = temperatureHandler->getTemperature(gridPosition,
					realTime);

			//xi = 4; ///FIXME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
				psiCluster = std::dynamic_pointer_cast < PSICluster
						> (network->get("He", i));
				computePartialsForDiffusion(psiCluster, temperature, sx, val,
						row, col, xi, xs, size);
				ierr = MatSetValuesLocal(*J, 1, row, 3, col, val, ADD_VALUES);
				checkPetscError(ierr);
			}

			/* 1V and 1I are the only other clusters that diffuse */

			// ----- Vacancy Diffusion -----
			psiCluster = std::dynamic_pointer_cast < PSICluster
					> (network->get("V", 1));
			if (psiCluster) {
				computePartialsForDiffusion(psiCluster, temperature, sx, val,
						row, col, xi, xs, size);
				ierr = MatSetValuesLocal(*J, 1, row, 3, col, val, ADD_VALUES);
				checkPetscError(ierr);
			}

			// ----- Interstitial Diffusion -----
			psiCluster = std::dynamic_pointer_cast < PSICluster
					> (network->get("I", 1));
			if (psiCluster) {
				computePartialsForDiffusion(psiCluster, temperature, sx, val,
						row, col, xi, xs, size);
				ierr = MatSetValuesLocal(*J, 1, row, 3, col, val, ADD_VALUES);
				checkPetscError(ierr);
			}
			//break;   // Uncomment this line for debugging in a single cell.
		}
		computeJacobianDiffusionTerms->stop();

		ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);
		checkPetscError(ierr);
		ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);
		checkPetscError(ierr);
//		ierr = MatSetOption(*J, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
//		checkPetscError(ierr);
		ierr = MatStoreValues(*J);
		checkPetscError(ierr);
//		MatSetFromOptions(*J);
		initialized = PETSC_TRUE;
		// Debug line for viewing the matrix
		//MatView(*J, PETSC_VIEWER_STDOUT_WORLD);
	} else {
		ierr = MatRetrieveValues(*J);
		checkPetscError(ierr);
	}

	/* ----- Compute the partial derivatives for the reaction term at each
	 * grid point ----- */

	// Create a new row array of size n
	PetscInt localPDColIds[size];
	PetscInt rowId = 0;
	// Create arrays for storing the partial derivatives
	std::vector<double> reactingPartialsForCluster(size, 0.0);
	std::vector<double> allPartialsForCluster;
	int pdColIdsVectorSize = 0;

	computeReactionTermPartials->start();
	// Loop over the grid points
	for (xi = xs; xi < xs + xm; xi++) {

		// Vector representing the position at which the flux will be calculated
		// Currently we are only in 1D
		std::vector<double> gridPosition = { 0, xi * hx, 0 };
		// Get the temperature handler that will be used to compute fluxes.
		auto temperatureHandler = PetscSolver::getTemperatureHandler();
		auto temperature = temperatureHandler->getTemperature(gridPosition,
				realTime);

		//xi = 4; ///FIXME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		// Copy data into the PSIClusterReactionNetwork so that it can
		// compute the new concentrations.
		concOffset = concs + size * xi;
		network->updateConcentrationsFromArray(concOffset);
		// Get the reactants
		reactants = network->getAll();
		updateJacobianCol->start();
		// Update the column in the Jacobian that represents each reactant
		for (int i = 0; i < size; i++) {
			reactant = reactants->at(i);
			// Get the reactant index
			reactantIndex = reactant->getId() - 1;
			// Get the column id
			rowId = (xi - xs + 1) * size + reactantIndex;
			// Get the partial derivatives
			allPartialsForCluster = reactant->getPartialDerivatives(
					temperature);
			// Set the row indices
			psiCluster = std::dynamic_pointer_cast < PSICluster > (reactant);
//			std::cout << xi << " " << xs << " " << size << " " << (xi - xs + 1)*size << std::endl;
//			std::cout << "PD for " << psiCluster->getName() << "_" << psiCluster->getSize() << " at " << reactantIndex << std::endl;
//			for (int k = 0; k < allPartialsForCluster.size(); k++) {
//				std::cout << "pd[" << k << "] = " << allPartialsForCluster[k] << std::endl;
//			}
			// Get the list of column ids from the map
			auto pdColIdsVector = dFillMap.at(reactantIndex);
			pdColIdsVectorSize = pdColIdsVector.size(); //Number of partial derivatives
			// Loop over the list of column ids
			for (int j = 0; j < pdColIdsVectorSize; j++) {
				// Calculate the appropriate index to match the dfill array configuration
				localPDColIds[j] = (xi - xs + 1) * size + pdColIdsVector[j];
				// Get the partial derivative from the array of all of the partials
				reactingPartialsForCluster[j] =
						allPartialsForCluster[pdColIdsVector[j]];
//				std::cout << "dp[" << j << "] = " << pdColIdsVector[j] << " , [r,c] = "<< "[" << rowId << "," << localPDColIds[j] << "] = " << reactingPartialsForCluster[j]<< std::endl;
			}
			// Update the matrix
			ierr = MatSetValuesLocal(*J, 1, &rowId, pdColIdsVectorSize,
					localPDColIds, reactingPartialsForCluster.data(),
					ADD_VALUES);
			checkPetscError(ierr);
		}
		updateJacobianCol->stop();
		// Uncomment this line for debugging in a single cell.
		//break;
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
	*str = SAME_NONZERO_PATTERN;
	ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);
	checkPetscError(ierr);
	ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);
	checkPetscError(ierr);

	// Boundary conditions
	// Loop over the grid points
	for (xi = xs; xi < xs + xm; xi++) {
		if (xi == 0) {
			// Get the reactants
			reactants = network->getAll();
			// Loop on the reactants
			for (int i = 0; i < size; i++) {
				reactant = reactants->at(i);
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
					if (pdColIdsVector[j] == reactantIndex)
						reactingPartialsForCluster[j] = 1.0;
					else
						reactingPartialsForCluster[j] = 0.0;
				}

				// Update the matrix
				ierr = MatSetValuesLocal(*J, 1, &rowId, pdColIdsVectorSize,
						localPDColIds, reactingPartialsForCluster.data(),
						INSERT_VALUES);
				checkPetscError(ierr);
			}
		}
	}

	// Assemble again
	ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);
	checkPetscError(ierr);
	ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);
	checkPetscError(ierr);

	if (*A != *J) {
		ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
		checkPetscError(ierr);
		ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
		checkPetscError(ierr);
	}
	PetscFunctionReturn(0);

}

PetscErrorCode callRHSJacobian(TS ts, PetscReal ftime, Vec C, Mat *A, Mat *J,
		MatStructure *str, void *ptr) {
	PetscErrorCode ierr;
	RHSJacobianTimer->start();
	ierr = RHSJacobian(ts, ftime, C, A, J, str, &ptr);
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
	std::shared_ptr<Reactant> reactant;

	// Fill the diagonal block if the sizes match up
	if (diagFillSize == size) {
		auto reactants = network->getAll();
		// Get the connectivity for each reactant
		for (i = 0; i < numReactants; i++) {
			// Get the reactant and its connectivity
			reactant = reactants->at(i);
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
	computeNewFluxes = handlerRegistry->getTimer("computeNewFluxes");
	computeIncidentFlux = handlerRegistry->getTimer("computeIncidentFlux");
	computeAllDiffusion = handlerRegistry->getTimer("computeAllDiffusion");
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

	// Debug
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

	// Check the network before getting busy.
	if (!network) {
		throw std::string("PetscSolver Exception: Network not set!");
	}

	// Set the output precision for std::out
	std::cout.precision(16);

	PetscFunctionBeginUser;
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = DMDACreate1d(PETSC_COMM_WORLD, DMDA_BOUNDARY_MIRROR, -8, dof, 1,
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
	std::shared_ptr<Reactant> reactant;
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
	ierr = TSSetInitialTimeStep(ts, 0.0, 1.0e-8);
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
