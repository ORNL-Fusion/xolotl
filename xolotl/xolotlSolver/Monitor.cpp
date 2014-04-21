// Includes
#include "PetscSolver.h"
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

namespace xolotlSolver {

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
 * This is a monitoring operation that displays He and V as a function of space
 * and cluster size for each time step. It is not a member variable of the class
 * because the monitoring code requires a C callback function (via a function
 * pointer).
 */
static PetscErrorCode monitorSolve(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int size = PetscSolver::getNetwork()->size();
	// The array of cluster names
	std::vector<std::string> names(size);
	// Header and output string streams
	std::stringstream header, outputData;
	std::ofstream outputFile;
	PetscErrorCode ierr;
	PetscReal *solutionArray, *gridPointSolution, x, hx;
	Vec localSolution;
	PetscInt xs, xm, Mx;
	int xi, i;

	PetscFunctionBeginUser;

	// Get the processor id
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Create a stream for naming the file
	std::stringstream outputFileNameStream;
	outputFileNameStream << "xolotl_out_" << procId << "_" << timestep;

	// Open the text output file
	outputFile.open(outputFileNameStream.str());

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the local vector, which is capital when running in parallel,
	// and put it into solutionArray
	ierr = DMGetLocalVector(da, &localSolution);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalBegin(da, solution, INSERT_VALUES, localSolution);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalEnd(da, solution, INSERT_VALUES, localSolution);
	checkPetscError(ierr);
	ierr = DMDAVecGetArray(da, localSolution, &solutionArray);
	checkPetscError(ierr);


	// Create the header for the file
	auto reactants = PetscSolver::getNetwork()->getAll();
	std::shared_ptr<PSICluster> cluster;
	header << "# t x ";
	for (int i = 0; i < size; i++) {
		// Get the cluster from the list, its id and composition
		cluster = std::dynamic_pointer_cast<PSICluster>(reactants->at(i));
		int id = cluster->getId() - 1;
		auto composition = cluster->getComposition();
		// Make the header entry
		std::stringstream name;
		name << (cluster->getName()).c_str() << "_(" << composition["He"] << ","
				<< composition["V"] << "," << composition["I"] << ") ";
		// Push the header entry on the array
		name >> names[id];
	}
	for (int i = 0; i < size; i++) {
		header << names[i] << " ";
	}
	header << "\n";

	// Writing the header in the output file
	outputFile << header.str();

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr);
	// Setup some step size variables
	hx = 8.0 / (PetscReal)(Mx - 1);
	checkPetscError(ierr);
	// Print the solution data
	for (xi = xs; xi < xs + xm; xi++) {
		// Dump x
		x = xi * hx;
		outputData << timestep << " " << x << " ";
		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray + size * xi;
		// Update the concentrations in the network to have physics results
		// (non negative)
		PetscSolver::getNetwork()->updateConcentrationsFromArray(
				gridPointSolution);
		// Get the concentrations from the network
		double concentrations[size];
		double * concentration = &concentrations[0];
		PetscSolver::getNetwork()->fillConcentrationsArray(concentration);
		// Dump the data to the stream
		for (i = 0; i < size; i++) {
			outputData << concentration[i] << " ";
		}
		// End the line
		outputData << "\n";
	}
	// Dump the data to file
	outputFile << outputData.str();
	// Restore the array and close the output file
	VecRestoreArray(solution, &solutionArray);
	outputFile.close();

	PetscFunctionReturn(0);
}

/**
 * This operation sets up a monitor that will display He as a function of space
 * and cluster size for each time step. It was not made a member function so
 * that it would be consistent with the other monitor callbacks.
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
PetscErrorCode setupPetscMonitor(TS ts) {
	PetscErrorCode ierr;
	PetscBool flg;

	PetscFunctionBeginUser;
	ierr = PetscOptionsHasName(NULL, "-mymonitor", &flg);
	checkPetscError(ierr);
	if (!flg)
		PetscFunctionReturn(0);

	ierr = TSMonitorSet(ts, monitorSolve, NULL, NULL);
	checkPetscError(ierr);
	PetscFunctionReturn(0);
}

}

/* end namespace xolotlSolver */
