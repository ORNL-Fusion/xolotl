// Includes
#include "PetscSolver.h"
#include <xolotlPerf.h>
#include <HDF5NetworkLoader.h>
#include <TemperatureHandler.h>
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
#include <DiffusionHandler.h>

using namespace xolotlCore;

/*
 C_t =  -D*C_xx + F(C) + R(C) + D(C) from Brian Wirth's SciDAC project.

 D*C_xx  - diffusion of He and V and I
 F(C)    - forcing function; He being created.
 R(C)    - reaction terms   (clusters combining)
 D(C)    - dissociation terms (cluster breaking up)

 Sample Options:
 -ts_monitor_draw_solution               -- plot the solution for each concentration as a function of x each in a separate 1d graph
 -draw_fields_by_name 1-He-2-V,1-He 	 -- only plot the solution for these two concentrations
 -da_grid_x <nx>						 -- number of grid points in the x direction
 -ts_max_steps <maxsteps>                -- maximum number of time-steps to take
 -ts_final_time <time>                   -- maximum time to compute to
 -ts_dt <size>							 -- initial size of the time step

 Rules for maximum number of He allowed for V in cluster
 */

namespace xolotlSolver {

//Timer for RHSFunction()
std::shared_ptr<xolotlPerf::ITimer> RHSFunctionTimer;

////Timer for RHSJacobian()
std::shared_ptr<xolotlPerf::ITimer> RHSJacobianTimer;

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
// Allocate the static diffusion handler
std::shared_ptr<IDiffusionHandler> PetscSolver::diffusionHandler;
// Allocate the static step size
double PetscSolver::hx;

extern PetscErrorCode RHSFunction(TS, PetscReal, Vec, Vec, void*);
extern PetscErrorCode RHSJacobian(TS, PetscReal, Vec, Mat, Mat);
extern PetscErrorCode setupPetscMonitor(TS);

TS ts; /* nonlinear solver */
Vec C; /* solution */
PetscErrorCode ierr;
DM da; /* manages the grid data */
PetscInt *ofill, *dfill;

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
 * A boolean that is true if the temperature has changed. It is set to true
 * in the RHSFunction, and back to false once the off-diagonal part of the
 * Jacobian is computed in the RHSJacobian method.
 */
static bool temperatureChanged = false;

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

	// Get the step size
	double hx = PetscSolver::getStepSize();

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
	std::shared_ptr<HDF5NetworkLoader> HDF5Loader = std::dynamic_pointer_cast
			< HDF5NetworkLoader > (networkLoader);
	auto fileName = HDF5Loader->getFilename();

	// Get the last time step written in the HDF5 file
	int tempTimeStep = -2;
	bool hasConcentrations = HDF5Utils::hasConcentrationGroup(fileName,
			tempTimeStep);

	// Loop on all the grid points
	for (i = xs; i < xs + xm; i++) {
		concOffset = concentrations + size * i;
		// Loop on all the clusters to initialize at 0.0
		for (int k = 0; k < size; k++) {
			concOffset[k] = 0.0;
		}

//		if (i > 0) {
//			int k = 14; // initial concentration for V only
//			concOffset[k] = 0.000315 / hx;
//		}

//		// Uncomment this for debugging
//		if (i > 0) {
//			for (int k = 0; k < size; k++) {
//				concOffset[k] = 0.001;
//			}
//		}
	}

	// If the concentration must be set from the HDF5 file
	if (hasConcentrations) {
		// Loop on the full grid
		for (int i = 0; i < Mx; i++) {
			// Read the concentrations from the HDF5 file
			auto concVector = HDF5Utils::readGridPoint(fileName, tempTimeStep,
					i);

			// Change the concentration only if we are on the locally owned part of the grid
			if (i >= xs && i < xs + xm) {
				concOffset = concentrations + size * i;
				// Loop on the concVector size
				for (int k = 0; k < concVector.size(); k++) {
					concOffset[(int) concVector.at(k).at(0)] =
							concVector.at(k).at(1);
				}
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
	// Start the RHSFunction Timer
	RHSFunctionTimer->start();

	// Important petsc stuff (related to the grid mostly)
	DM da;
	PetscErrorCode ierr;
	PetscInt xi, Mx, xs, xm;
	PetscReal x;
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
	// Get the temperature handler that will be used to compute fluxes.
	auto temperatureHandler = PetscSolver::getTemperatureHandler();
	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = PetscSolver::getFluxHandler();
	auto incidentFluxVector = fluxHandler->getIncidentFluxVec();
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

	// Setup some step size variables
	double hx = PetscSolver::getStepSize();
	double sx = 1.0 / (hx * hx);

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

	// Get the diffusion handler
	auto diffusionHandler = PetscSolver::getDiffusionHandler();

	// Loop over grid points computing ODE terms for each grid point
	size = network->size();
	for (xi = xs; xi < xs + xm; xi++) {

//		xi = 1; // Uncomment this line for debugging in a single cell.

		x = xi * hx;

		// Vector representing the position at which the flux will be calculated
		// Currently we are only in 1D
		std::vector<double> gridPosition = { x, 0, 0 };
		auto temperature = temperatureHandler->getTemperature(gridPosition,
				realTime);

		// Update the network if the temperature changed
		if (!xolotlCore::equal(temperature, lastTemperature)) {
			network->setTemperature(temperature);
			lastTemperature = temperature;
			temperatureChanged = true;
		}

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
		heCluster = (PSICluster *) network->get("He", 1);
		if (heCluster) {
			reactantIndex = heCluster->getId() - 1;
			// Update the concentration of the cluster
			updatedConcOffset[reactantIndex] += incidentFluxVector[xi];
		}

		// ---- Compute diffusion over the locally owned part of the grid -----
		diffusionHandler->computeDiffusion(network, sx, concOffset,
				leftConcOffset, rightConcOffset, updatedConcOffset);

		// ----- Compute all of the new fluxes -----
		for (int i = 0; i < size; i++) {
			cluster = (PSICluster *) allReactants->at(i);
			// Compute the flux
			flux = cluster->getTotalFlux();
			// Update the concentration of the cluster
			reactantIndex = cluster->getId() - 1;
			updatedConcOffset[reactantIndex] += flux;
//			std::cout << "New flux = " << flux << " "
//					<< cluster->getConcentration() << std::endl;
		}

		// Boundary conditions
		if (xi == 0) {
			for (int i = 0; i < size; i++) {
				updatedConcOffset[i] = 1.0 * concs[i];
			}
		}

		// Uncomment this line for debugging in a single cell.
//		break;
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArray(da, localC, &concs);
	checkPetscError(ierr);
	ierr = DMDAVecRestoreArray(da, F, &updatedConcs);
	checkPetscError(ierr);
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr);

	// Stop the RHSFunction Timer
	RHSFunctionTimer->stop();

	PetscFunctionReturn(0);

}

PetscErrorCode callRHSFunction(TS ts, PetscReal ftime, Vec C, Vec F,
		void *ptr) {

	PetscErrorCode ierr;
	ierr = RHSFunction(ts, ftime, C, F, &ptr);

	return ierr;
}

#undef __FUNCT__
#define __FUNCT__ "RHSJacobian"
/*
 Compute the Jacobian entries based on IFuction() and insert them into the matrix
 */
PetscErrorCode RHSJacobian(TS ts, PetscReal ftime, Vec C, Mat A, Mat J,
		void *ptr) {
	// Start the RHSJacobian timer
	RHSJacobianTimer->start();

	DM da;
	PetscErrorCode ierr;
	PetscInt xi, Mx, xs, xm, i;
	PetscReal *concs, *updatedConcs;
	double * concOffset;
	Vec localC;
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

	// Setup some step size variables
	double hx = PetscSolver::getStepSize();
	double sx = 1.0 / (hx * hx);

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
	checkPetscError(ierr);

	// Only compute the linear part of the Jacobian if the temperature has changed
	if (temperatureChanged) {

		// Create the diffusion handler
		auto diffusionHandler = PetscSolver::getDiffusionHandler();
		// Get the total number of diffusing clusters
		const int nDiff = diffusionHandler->getNumberOfDiffusing();

		/*
		 Loop over grid points computing Jacobian terms for diffusion at each
		 grid point
		 */
		for (xi = xs; xi < xs + xm; xi++) {

//			xi = 1; // Uncomment this line for debugging in a single cell

			// Boundary conditions
			if (xi == 0) continue;

			// Copy data into the PSIClusterReactionNetwork so that it can
			// compute the new concentrations.
			concOffset = concs + size * xi;
			network->updateConcentrationsFromArray(concOffset);

			// Initialize the rows, columns, and values to set in the Jacobian
			PetscInt row[nDiff], col[3*nDiff];
			PetscReal val[3*nDiff];
			// Get the pointer on them for the compute diffusion method
			PetscInt *rowPointer = &row[0];
			PetscInt *colPointer = &col[0];
			PetscReal *valPointer = &val[0];

			/* -------------------------------------------------------------
			 ---- Compute diffusion over the locally owned part of the grid
			 */
			diffusionHandler->computePartialsForDiffusion(network, sx, valPointer,
					rowPointer, colPointer, xi, xs);

			// Loop on the number of diffusion cluster to set the values in the Jacobian
			for(int i = 0; i < nDiff; i++) {
				rowPointer = &row[i];
				colPointer = &col[3*i];
				valPointer = &val[3*i];

				ierr = MatSetValuesLocal(J, 1, rowPointer, 3, colPointer, valPointer, ADD_VALUES);
				checkPetscError(ierr);
			}

//			break;   // Uncomment this line for debugging in a single cell.
		}

		ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
		checkPetscError(ierr);
		ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
		checkPetscError(ierr);
//		ierr = MatSetOption(J, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
//		checkPetscError(ierr);
		ierr = MatStoreValues(J);
		checkPetscError(ierr);
//		MatSetFromOptions(J);
		temperatureChanged = false;
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

	// Loop over the grid points
	for (xi = xs; xi < xs + xm; xi++) {

//		xi = 1; // Uncomment this line for debugging in a single cell

		// Copy data into the PSIClusterReactionNetwork so that it can
		// compute the new concentrations.
		concOffset = concs + size * xi;
		network->updateConcentrationsFromArray(concOffset);
		// Update the column in the Jacobian that represents each reactant
		for (int i = 0; i < size; i++) {
			auto reactant = allReactants->at(i);
			// Get the reactant index
			reactantIndex = reactant->getId() - 1;
			// Get the column id
			rowId = (xi - xs + 1) * size + reactantIndex;
			// Get the partial derivatives
			reactant->getPartialDerivatives(clusterPartials);
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
		// Uncomment this line for debugging in a single cell.
//		break;
	}

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
	if (xs == 0) {
		// Loop on the reactants
		for (int i = 0; i < size; i++) {
			auto reactant = allReactants->at(i);
			// Get the reactant index
			reactantIndex = reactant->getId() - 1;
			// Get the row id
			rowId = size + reactantIndex;

			// Get the list of column ids from the map
			auto pdColIdsVector = dFillMap.at(reactantIndex);
			pdColIdsVectorSize = pdColIdsVector.size(); //Number of partial derivatives
			// Loop over the list of column ids
			for (int j = 0; j < pdColIdsVectorSize; j++) {
				// Calculate the appropriate index to match the dfill array configuration
				localPDColIds[j] = size + pdColIdsVector[j];
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

	// Stop the RHSJacobian timer
	RHSJacobianTimer->stop();

	PetscFunctionReturn(0);

}

PetscErrorCode callRHSJacobian(TS ts, PetscReal ftime, Vec C, Mat A, Mat J,
		void *ptr) {

	PetscErrorCode ierr;
	ierr = RHSJacobian(ts, ftime, C, A, J, &ptr);

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

	RHSFunctionTimer = handlerRegistry->getTimer("RHSFunctionTimer");
	RHSJacobianTimer = handlerRegistry->getTimer("RHSJacobianTimer");
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
	if (procId == 0) {
		std::cout << "\nPETScSolver Message: " << "Master loaded network of size "
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
		std::shared_ptr<ITemperatureHandler> temperatureHandler,
		double stepSize) {

	// Set the flux handler
	PetscSolver::fluxHandler = fluxHandler;

	// Set the temperature handler
	PetscSolver::temperatureHandler = temperatureHandler;

	// Set the diffusion handler
	PetscSolver::diffusionHandler = std::make_shared<DiffusionHandler>();

	// Set the grid step size
	PetscSolver::hx = stepSize;
	
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
	std::shared_ptr<HDF5NetworkLoader> HDF5Loader = std::dynamic_pointer_cast
			< HDF5NetworkLoader > (networkLoader);
	auto fileName = HDF5Loader->getFilename();

	// Get starting conditions from HDF5 file
	int gridLength = 0;
	double time = 0.0, deltaTime = 1.0e-12;
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

	/*  The only spatial coupling in the Jacobian is due to diffusion.
	 *  The ofill (thought of as a dof by dof 2d (row-oriented) array represents
	 *  the nonzero coupling between degrees of freedom at one point with degrees
	 *  of freedom on the adjacent point to the left or right. A 1 at i,j in the
	 *  ofill array indicates that the degree of freedom i at a point is coupled
	 *  to degree of freedom j at the adjacent point.
	 *  In this case ofill has only a few diagonal entries since the only spatial
	 *  coupling is regular diffusion.
	 */
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &ofill);
	checkPetscError(ierr);
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &dfill);
	checkPetscError(ierr);
	ierr = PetscMemzero(ofill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr);
	ierr = PetscMemzero(dfill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr);

	// Fill ofill, the matrix of "off-diagonal" elements that represents diffusion
	PetscSolver::diffusionHandler->initializeOFill(network, ofill);

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
	checkPetscError(ierr);
	ierr = TSSetRHSJacobian(ts, NULL, NULL, callRHSJacobian, NULL);
	checkPetscError(ierr);
	ierr = TSSetSolution(ts, C);
	checkPetscError(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set solver options
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = TSSetInitialTimeStep(ts, time, deltaTime);
	checkPetscError(ierr);
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
	if (ts != NULL && C != NULL) {
		ierr = TSSolve(ts, C);
		checkPetscError(ierr);
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
