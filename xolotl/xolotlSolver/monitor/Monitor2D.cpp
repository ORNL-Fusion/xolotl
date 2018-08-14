// Includes
#include "PetscSolver.h"
#include <xolotlPerf.h>
#include <VizHandlerRegistryFactory.h>
#include <PlotType.h>
#include <CvsXDataProvider.h>
#include <CvsXYDataProvider.h>
#include <LabelProvider.h>
#include <Constants.h>
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <PSISuperCluster.h>
#include <NESuperCluster.h>
#include <MathUtils.h>
#include "RandomNumberGenerator.h"
#include "xolotlCore/io/XFile.h"
#include "xolotlSolver/monitor/Monitor.h"

namespace xolotlSolver {

// Declaration of the functions defined in Monitor.cpp
extern PetscErrorCode checkTimeStep(TS ts);
extern PetscErrorCode monitorTime(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode computeFluence(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode monitorPerf(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);

// Declaration of the variables defined in Monitor.cpp
extern std::shared_ptr<xolotlViz::IPlot> perfPlot;
extern double previousTime;
extern double timeStepThreshold;

//! How often HDF5 file is written
PetscReal hdf5Stride2D = 0.0;
//! Previous time for HDF5
PetscInt hdf5Previous2D = 0;
//! HDF5 output file name
std::string hdf5OutputName2D = "xolotlStop.h5";
//! The pointer to the 2D plot used in MonitorSurface.
std::shared_ptr<xolotlViz::IPlot> surfacePlot2D;
//! The variable to store the interstitial flux at the previous time step.
std::vector<double> previousIFlux2D;
//! The variable to store the total number of interstitials going through the surface.
std::vector<double> nInterstitial2D;
//! The variable to store the helium flux at the previous time step.
std::vector<double> previousHeFlux2D;
//! The variable to store the total number of helium going through the bottom.
std::vector<double> nHelium2D;
//! The variable to store the deuterium flux at the previous time step.
std::vector<double> previousDFlux2D;
//! The variable to store the total number of deuterium going through the bottom.
std::vector<double> nDeuterium2D;
//! The variable to store the tritium flux at the previous time step.
std::vector<double> previousTFlux2D;
//! The variable to store the total number of tritium going through the bottom.
std::vector<double> nTritium2D;
//! The variable to store the sputtering yield at the surface.
double sputteringYield2D = 0.0;
// The vector of depths at which bursting happens
std::vector<std::pair<int, int> > depthPositions2D;

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "startStop2D")
/**
 * This is a monitoring method that will update an hdf5 file at each time step.
 */
PetscErrorCode startStop2D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *) {
	// Initial declaration
	PetscErrorCode ierr;
	const double ***solutionArray, *gridPointSolution;
	PetscInt xs, xm, Mx, ys, ym, My;

	PetscFunctionBeginUser;

	// Compute the dt
	double dt = time - previousTime;

	// Don't do anything if it is not on the stride
	if ((int) ((time + dt / 10.0) / hdf5Stride2D) <= hdf5Previous2D)
		PetscFunctionReturn(0);

	// Update the previous time
	hdf5Previous2D++;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto& network = solverHandler.getNetwork();

	// Network size
	const int dof = network.getDOF();

	// Create an array for the concentration
	double concArray[dof][2];

	// Get the vector of positions of the surface
	std::vector<int> surfaceIndices;
	for (PetscInt i = 0; i < My; i++) {
		surfaceIndices.push_back(solverHandler.getSurfacePosition(i));
	}

	// Open the existing checkpoint file.
	xolotlCore::XFile checkpointFile(hdf5OutputName2D, PETSC_COMM_WORLD,
			xolotlCore::XFile::AccessMode::OpenReadWrite);

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	CHKERRQ(ierr);

	// Add a concentration sub group
	auto concGroup = checkpointFile.getGroup<
			xolotlCore::XFile::ConcentrationGroup>();
	assert(concGroup);
	auto tsGroup = concGroup->addTimestepGroup(timestep, time, previousTime,
			currentTimeStep);

	// Write the surface positions and the associated interstitial quantities
	// in the concentration sub group
	tsGroup->writeSurface2D(surfaceIndices, nInterstitial2D, previousIFlux2D);

	// Write the bottom impurity information if the bottom is a free surface
	if (solverHandler.getRightOffset() == 1)
		tsGroup->writeBottom2D(nHelium2D, previousHeFlux2D, nDeuterium2D,
				previousDFlux2D, nTritium2D, previousTFlux2D);

	// Loop on the full grid
	for (PetscInt j = 0; j < My; j++) {
		for (PetscInt i = 0; i < Mx; i++) {
			// Size of the concentration that will be stored
			int concSize = -1;
			// To save which proc has the information
			int concId = 0;

			// If it is the locally owned part of the grid
			if (i >= xs && i < xs + xm && j >= ys && j < ys + ym) {
				// Get the pointer to the beginning of the solution data for this grid point
				gridPointSolution = solutionArray[j][i];

				// Loop on the concentrations
				for (int l = 0; l < dof; l++) {
					if (gridPointSolution[l] > 1.0e-16
							|| gridPointSolution[l] < -1.0e-16) {
						// Increase concSize
						concSize++;
						// Fill the concArray
						concArray[concSize][0] = (double) l;
						concArray[concSize][1] = gridPointSolution[l];
					}
				}

				// Increase concSize one last time
				concSize++;

				// Save the procId
				concId = procId;
			}

			// Get which processor will send the information
			int concProc = 0;
			MPI_Allreduce(&concId, &concProc, 1, MPI_INT, MPI_SUM,
					PETSC_COMM_WORLD);

			// Broadcast the size
			MPI_Bcast(&concSize, 1, MPI_INT, concProc, PETSC_COMM_WORLD);

			// Skip the grid point if the size is 0
			if (concSize == 0)
				continue;

			// Transfer the data everywhere from the local grid
			MPI_Bcast(&(concArray[0][0]), 2 * concSize, MPI_DOUBLE, concProc,
					PETSC_COMM_WORLD);

			// All processes create the dataset and fill it
			tsGroup->writeConcentrationDataset(concSize, concArray, i, j);
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeHeliumRetention2D")
/**
 * This is a monitoring method that will compute the helium retention
 */
PetscErrorCode computeHeliumRetention2D(TS ts, PetscInt, PetscReal time,
		Vec solution, void *) {
	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm, ys, ym, Mx, My;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = solverHandler.getFluxHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid in the x direction
	auto grid = solverHandler.getXGrid();

	// Setup step size variables
	double hy = solverHandler.getStepSizeY();

	// Get the network
	auto& network = solverHandler.getNetwork();

	// Get the array of concentration
	double ***solutionArray, *gridPointSolution;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration over the grid
	double heConcentration = 0.0, dConcentration = 0.0, tConcentration = 0.0;

	// Loop on the grid
	for (PetscInt yj = ys; yj < ys + ym; yj++) {
		// Get the surface position
		int surfacePos = solverHandler.getSurfacePosition(yj);

		for (PetscInt xi = xs; xi < xs + xm; xi++) {

			// Boundary conditions
			if (xi < surfacePos || xi == Mx - 1)
				continue;

			// Get the pointer to the beginning of the solution data for this grid point
			gridPointSolution = solutionArray[yj][xi];

			// Update the concentration in the network
			network.updateConcentrationsFromArray(gridPointSolution);

			// Get the total atom concentrations at this grid point
			heConcentration += network.getTotalAtomConcentration(0)
					* (grid[xi + 1] - grid[xi]) * hy;
			dConcentration += network.getTotalAtomConcentration(1)
					* (grid[xi + 1] - grid[xi]) * hy;
			tConcentration += network.getTotalAtomConcentration(2)
					* (grid[xi + 1] - grid[xi]) * hy;
		}
	}

	// Get the current process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Sum all the concentrations through MPI reduce
	double totalHeConcentration = 0.0;
	MPI_Reduce(&heConcentration, &totalHeConcentration, 1, MPI_DOUBLE, MPI_SUM,
			0, PETSC_COMM_WORLD);
	double totalDConcentration = 0.0;
	MPI_Reduce(&dConcentration, &totalDConcentration, 1, MPI_DOUBLE, MPI_SUM, 0,
			PETSC_COMM_WORLD);
	double totalTConcentration = 0.0;
	MPI_Reduce(&tConcentration, &totalTConcentration, 1, MPI_DOUBLE, MPI_SUM, 0,
			PETSC_COMM_WORLD);

	// Get the total size of the grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	CHKERRQ(ierr);

	// Look at the fluxes going in the bulk if the bottom is a free surface
	if (solverHandler.getRightOffset() == 1) {
		// Set the bottom surface position
		int xi = Mx - 2;

		// Loop on every Y position
		for (PetscInt j = 0; j < My; j++) {

			// Value to know on which processor is the bottom
			int bottomProc = 0;

			// Check we are on the right proc
			if (xi >= xs && xi < xs + xm && j >= ys && j < ys + ym) {
				// Get the delta time from the previous timestep to this timestep
				double dt = time - previousTime;
				// Compute the total number of impurities that went in the bulk
				nHelium2D[j] += previousHeFlux2D[j] * dt;
				nDeuterium2D[j] += previousDFlux2D[j] * dt;
				nTritium2D[j] += previousTFlux2D[j] * dt;

				// Get the pointer to the beginning of the solution data for this grid point
				gridPointSolution = solutionArray[j][xi];

				// Factor for finite difference
				double hxLeft = grid[xi + 1] - grid[xi];
				double hxRight = grid[xi + 2] - grid[xi + 1];
				double factor = 2.0 / (hxLeft + hxRight);

				// Initialize the value for the flux
				double newFlux = 0.0;
				// Consider each helium cluster.
				for (auto const& heMapItem : network.getAll(ReactantType::He)) {
					// Get the cluster
					auto const& cluster = *(heMapItem.second);
					// Get its id and concentration
					int id = cluster.getId() - 1;
					double conc = gridPointSolution[id];
					// Get its size and diffusion coefficient
					int size = cluster.getSize();
					double coef = cluster.getDiffusionCoefficient(xi);
					// Compute the flux going to the right
					newFlux += (double) size * factor * coef * conc;
				}
				// Update the helium flux
				previousHeFlux2D[j] = newFlux;

				// Initialize the value for the flux
				newFlux = 0.0;
				// Consider each deuterium cluster.
				for (auto const& dMapItem : network.getAll(ReactantType::D)) {
					// Get the cluster
					auto const& cluster = *(dMapItem.second);
					// Get its id and concentration
					int id = cluster.getId() - 1;
					double conc = gridPointSolution[id];
					// Get its size and diffusion coefficient
					int size = cluster.getSize();
					double coef = cluster.getDiffusionCoefficient(xi);
					// Compute the flux going to the right
					newFlux += (double) size * factor * coef * conc;
				}
				// Update the deuterium flux
				previousDFlux2D[j] = newFlux;

				// Initialize the value for the flux
				newFlux = 0.0;
				// Consider each tritium cluster.
				for (auto const& tMapItem : network.getAll(ReactantType::T)) {
					// Get the cluster
					auto const& cluster = *(tMapItem.second);
					// Get its id and concentration
					int id = cluster.getId() - 1;
					double conc = gridPointSolution[id];
					// Get its size and diffusion coefficient
					int size = cluster.getSize();
					double coef = cluster.getDiffusionCoefficient(xi);
					// Compute the flux going to the right
					newFlux += (double) size * factor * coef * conc;
				}
				// Update the tritium flux
				previousTFlux2D[j] = newFlux;

				// Set the bottom processor
				bottomProc = procId;
			}

			// Get which processor will send the information
			int bottomId = 0;
			MPI_Allreduce(&bottomProc, &bottomId, 1, MPI_INT, MPI_SUM,
					PETSC_COMM_WORLD);

			// Send the information about impurities
			// to the other processes
			MPI_Bcast(&nHelium2D[j], 1, MPI_DOUBLE, bottomId, PETSC_COMM_WORLD);
			MPI_Bcast(&previousHeFlux2D[j], 1, MPI_DOUBLE, bottomId,
					PETSC_COMM_WORLD);
			MPI_Bcast(&nDeuterium2D[j], 1, MPI_DOUBLE, bottomId,
					PETSC_COMM_WORLD);
			MPI_Bcast(&previousDFlux2D[j], 1, MPI_DOUBLE, bottomId,
					PETSC_COMM_WORLD);
			MPI_Bcast(&nTritium2D[j], 1, MPI_DOUBLE, bottomId,
					PETSC_COMM_WORLD);
			MPI_Bcast(&previousTFlux2D[j], 1, MPI_DOUBLE, bottomId,
					PETSC_COMM_WORLD);
		}
	}

	// Master process
	if (procId == 0) {
		// Compute the total surface irradiated by the helium flux
		double surface = (double) My * hy;

		// Rescale the concentration
		totalHeConcentration = totalHeConcentration / surface;
		totalDConcentration = totalDConcentration / surface;
		totalTConcentration = totalTConcentration / surface;
		double totalHeBulk = 0.0, totalDBulk = 0.0, totalTBulk = 0.0;
		// Look if the bottom is a free surface
		if (solverHandler.getRightOffset() == 1) {
			for (int i = 0; i < My; i++) {
				totalHeBulk += nHelium2D[i];
				totalDBulk += nDeuterium2D[i];
				totalTBulk += nTritium2D[i];
			}
			totalHeBulk = totalHeBulk / surface;
			totalDBulk = totalDBulk / surface;
			totalTBulk = totalTBulk / surface;
		}

		// Get the fluence
		double fluence = fluxHandler->getFluence();

		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		std::cout << "Helium content = " << totalHeConcentration << std::endl;
		std::cout << "Deuterium content = " << totalDConcentration << std::endl;
		std::cout << "Tritium content = " << totalTConcentration << std::endl;
		std::cout << "Fluence = " << fluence << "\n" << std::endl;

		// Uncomment to write the retention and the fluence in a file
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt", ios::app);
		outputFile << fluence << " " << totalHeConcentration << " "
				<< totalDConcentration << " " << totalTConcentration << " "
				<< totalHeBulk << " " << totalDBulk << " " << totalTBulk
				<< std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeTRIDYN2D")
/**
 * This is a monitoring method that will compute the data to send to TRIDYN
 */
PetscErrorCode computeTRIDYN2D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm, ys, ym;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto& network = solverHandler.getNetwork();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	CHKERRQ(ierr);

	// Get the total size of the grid rescale the concentrations
	PetscInt Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// Get the array of concentration
	double ***solutionArray, *gridPointSolution;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Create the output file
	std::ofstream outputFile;
	if (procId == 0) {
		std::stringstream name;
		name << "TRIDYN_" << timestep << ".dat";
		outputFile.open(name.str());
	}

	// Loop on the entire grid
	for (int xi = 0; xi < Mx; xi++) {
		// Set x
		double x = grid[xi + 1] - grid[1];

		// Initialize the concentrations at this grid point
		double heLocalConc = 0.0, dLocalConc = 0.0, tLocalConc = 0.0,
				vLocalConc = 0.0, iLocalConc = 0.0;

		// Loop on the y
		for (PetscInt yj = ys; yj < ys + ym; yj++) {
			// Get the surface position
			int surfacePos = solverHandler.getSurfacePosition(yj);
			// Boundary conditions
			if (xi < surfacePos)
				continue;

			// If it is the locally owned part of the grid
			if (xi >= xs && xi < xs + xm) {
				// Get the pointer to the beginning of the solution data for this grid point
				gridPointSolution = solutionArray[yj][xi];

				// Update the concentration in the network
				network.updateConcentrationsFromArray(gridPointSolution);

				// Get the total helium concentration at this grid point
				heLocalConc += network.getTotalAtomConcentration(0);
				dLocalConc += network.getTotalAtomConcentration(1);
				tLocalConc += network.getTotalAtomConcentration(2);
				vLocalConc += network.getTotalVConcentration();
				iLocalConc += network.getTotalIConcentration();
			}
		}

		double heConc = 0.0, dConc = 0.0, tConc = 0.0, vConc = 0.0, iConc = 0.0;
		MPI_Reduce(&heLocalConc, &heConc, 1, MPI_DOUBLE, MPI_SUM, 0,
				PETSC_COMM_WORLD);
		MPI_Reduce(&dLocalConc, &dConc, 1, MPI_DOUBLE, MPI_SUM, 0,
				PETSC_COMM_WORLD);
		MPI_Reduce(&tLocalConc, &tConc, 1, MPI_DOUBLE, MPI_SUM, 0,
				PETSC_COMM_WORLD);
		MPI_Reduce(&vLocalConc, &vConc, 1, MPI_DOUBLE, MPI_SUM, 0,
				PETSC_COMM_WORLD);
		MPI_Reduce(&iLocalConc, &iConc, 1, MPI_DOUBLE, MPI_SUM, 0,
				PETSC_COMM_WORLD);

		// The master process writes computes the cumulative value and writes in the file
		if (procId == 0) {
			outputFile
					<< x
							- (grid[solverHandler.getSurfacePosition(0) + 1]
									- grid[1]) << " " << heConc / My << " "
					<< dConc / My << " " << tConc / My << " " << vConc / My
					<< " " << iConc / My << std::endl;
		}
	}

	// Close the file
	if (procId == 0) {
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorSurface2D")
/**
 * This is a monitoring method that will save 2D plots of the concentration of
 * a specific cluster at each grid point.
 */
PetscErrorCode monitorSurface2D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *) {
	// Initial declarations
	PetscErrorCode ierr;
	const double ***solutionArray, *gridPointSolution;
	PetscInt xs, xm, Mx, ys, ym, My;
	double x = 0.0, y = 0.0;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % 10 != 0)
		PetscFunctionReturn(0);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto& network = solverHandler.getNetwork();

	// Get the physical grid in the x direction
	auto grid = solverHandler.getXGrid();

	// Setup step size variables
	double hy = solverHandler.getStepSizeY();

	// Choice of the cluster to be plotted
	int iCluster = 19;

	// Create a Point vector to store the data to give to the data provider
	// for the visualization
	auto myPoints = std::make_shared<std::vector<xolotlViz::Point> >();
	// Create a point here so that it is not created and deleted in the loop
	xolotlViz::Point thePoint;

	// Loop on the full grid
	for (PetscInt j = 0; j < My; j++) {
		for (PetscInt i = 0; i < Mx; i++) {
			// If it is the locally owned part of the grid
			if (i >= xs && i < xs + xm && j >= ys && j < ys + ym) {
				// Get the pointer to the beginning of the solution data for this grid point
				gridPointSolution = solutionArray[j][i];
				// Compute x and y
				x = grid[i + 1] - grid[1];
				y = (double) j * hy;

				// If it is procId 0 just store the value in the myPoints vector
				if (procId == 0) {
					// Modify the Point with the gridPointSolution[iCluster] as the value
					// and add it to myPoints
					thePoint.value = gridPointSolution[iCluster];
					thePoint.t = time;
					thePoint.x = x;
					thePoint.y = y;
					myPoints->push_back(thePoint);
				}
				// Else, the values must be sent to procId 0
				else {
					// Send the value of the local position to the master process
					MPI_Send(&x, 1, MPI_DOUBLE, 0, 10, PETSC_COMM_WORLD);
					// Send the value of the local position to the master process
					MPI_Send(&y, 1, MPI_DOUBLE, 0, 11, PETSC_COMM_WORLD);

					// Send the value of the concentration to the master process
					MPI_Send(&gridPointSolution[iCluster], 1, MPI_DOUBLE, 0, 12,
							PETSC_COMM_WORLD);
				}
			}
			// Else if it is NOT the locally owned part of the grid but still procId == 0,
			// it should receive the values for the point and add them to myPoint
			else if (procId == 0) {
				// Get the position
				MPI_Recv(&x, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 10,
						PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&y, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 11,
						PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

				// and the concentration
				double conc = 0.0;
				MPI_Recv(&conc, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 12,
						PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

				// Modify the Point with the received values and add it to myPoints
				thePoint.value = conc;
				thePoint.t = time;
				thePoint.x = x;
				thePoint.y = y;
				myPoints->push_back(thePoint);
			}

			// Wait for everybody at each grid point
			MPI_Barrier(PETSC_COMM_WORLD);
		}
	}

	// Plot everything from procId == 0
	if (procId == 0) {
		// Get the data provider and give it the points
		surfacePlot2D->getDataProvider()->setPoints(myPoints);

		// Get the iCluster cluster to have access to its name
		auto const& reactants = network.getAll();
		IReactant const& cluster = reactants.at(iCluster);

		// Change the title of the plot and the name of the data
		std::stringstream title;
		title << cluster.getName();
		surfacePlot2D->getDataProvider()->setDataName(title.str());
		title << " concentration";
		surfacePlot2D->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		surfacePlot2D->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		surfacePlot2D->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << cluster.getName() << "_surface_TS" << timestep << ".png";
		surfacePlot2D->write(fileName.str());
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "eventFunction2D")
/**
 * This is a method that checks if the surface should move or bursting happen
 */
PetscErrorCode eventFunction2D(TS ts, PetscReal time, Vec solution,
		PetscScalar *fvalue, void *) {
	// Initial declaration
	PetscErrorCode ierr;
	double ***solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi, Mx, ys, ym, yj, My;
	fvalue[0] = 1.0, fvalue[1] = 1.0;
	depthPositions2D.clear();

	PetscFunctionBeginUser;

	// Gets the process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	CHKERRQ(ierr);

	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto& network = solverHandler.getNetwork();

	// Get the physical grid
	auto grid = solverHandler.getXGrid();
	// Get the step size in Y
	double hy = solverHandler.getStepSizeY();

	// Get the flux handler to know the flux amplitude.
	auto fluxHandler = solverHandler.getFluxHandler();
	double heliumFluxAmplitude = fluxHandler->getFluxAmplitude();

	// Get the delta time from the previous timestep to this timestep
	double dt = time - previousTime;

	// Work of the moving surface first
	if (solverHandler.moveSurface()) {
		// Write the initial surface positions
		if (procId == 0 && xolotlCore::equal(time, 0.0)) {
			std::ofstream outputFile;
			outputFile.open("surface.txt", ios::app);
			outputFile << time << " ";

			// Loop on the possible yj
			for (yj = 0; yj < My; yj++) {
				// Get the position of the surface at yj
				int surfacePos = solverHandler.getSurfacePosition(yj);
				outputFile << grid[surfacePos + 1] - grid[1] << " ";
			}
			outputFile << std::endl;
			outputFile.close();
		}

		// Get the initial vacancy concentration
		double initialVConc = solverHandler.getInitialVConc();

		// Loop on the possible yj
		for (yj = 0; yj < My; yj++) {

			// Compute the total density of intersitials that escaped from the
			// surface since last timestep using the stored flux
			nInterstitial2D[yj] += previousIFlux2D[yj] * dt;

			// Remove the sputtering yield since last timestep
			nInterstitial2D[yj] -= sputteringYield2D * heliumFluxAmplitude * dt
					* hy;

			// Get the position of the surface at yj
			const int surfacePos = solverHandler.getSurfacePosition(yj);
			xi = surfacePos + 1;
			const double hxLeft = grid[xi + 1] - grid[xi];
			const double hxRight = grid[xi + 2] - grid[xi + 1];

			// Initialize the value for the flux
			double newFlux = 0.0;

			// if xi is on this process
			if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
				// Get the concentrations at xi = surfacePos + 1
				gridPointSolution = solutionArray[yj][xi];

				// Factor for finite difference
				double factor = 2.0 / (hxLeft + hxRight);

				// Loop on all the interstitial clusters to add the contribution from deeper
				for (auto const& iMapItem : network.getAll(ReactantType::I)) {
					// Get the cluster
					auto const& cluster = *(iMapItem.second);
					// Get its id and concentration
					int id = cluster.getId() - 1;
					double conc = gridPointSolution[id];
					// Get its size and diffusion coefficient
					int size = cluster.getSize();
					double coef = cluster.getDiffusionCoefficient(xi);
					// Compute the flux going to the left
					newFlux += (double) size * factor * coef * conc * hy;
				}
			}

			// Check if the surface on the left and/or right sides are higher than at this grid point
			int yLeft = yj - 1, yRight = yj + 1;
			if (yLeft < 0)
				yLeft = My - 1; // Periodicity
			if (yRight == My)
				yRight = 0; // Periodicity
			// We want the position just at the surface now
			xi = surfacePos;
			// Do the left side first
			if (solverHandler.getSurfacePosition(yLeft) > surfacePos) {
				// if we are on the right process
				if (xi >= xs && xi < xs + xm && yLeft >= ys
						&& yLeft < ys + ym) {
					// Get the concentrations at xi = surfacePos
					gridPointSolution = solutionArray[yLeft][xi];

					// Loop on all the interstitial clusters to add the contribution from the left side
					for (auto const& iMapItem : network.getAll(ReactantType::I)) {
						// Get the cluster
						auto const& cluster = *(iMapItem.second);
						// Get its id and concentration
						int id = cluster.getId() - 1;
						double conc = gridPointSolution[id];
						// Get its size and diffusion coefficient
						int size = cluster.getSize();
						double coef = cluster.getDiffusionCoefficient(xi);
						// Compute the flux
						newFlux += ((double) size * coef * conc * hxLeft) / hy;
					}
				}
			}
			// Now do the right side
			if (solverHandler.getSurfacePosition(yRight) > surfacePos) {
				// if we are on the right process
				if (xi >= xs && xi < xs + xm && yRight >= ys
						&& yRight < ys + ym) {
					// Get the concentrations at xi = surfacePos + 1
					gridPointSolution = solutionArray[yRight][xi];

					// Loop on all the interstitial clusters to add the contribution from the left side
					for (auto const& iMapItem : network.getAll(ReactantType::I)) {
						// Get the cluster
						auto const& cluster = *(iMapItem.second);
						// Get its id and concentration
						int id = cluster.getId() - 1;
						double conc = gridPointSolution[id];
						// Get its size and diffusion coefficient
						int size = cluster.getSize();
						double coef = cluster.getDiffusionCoefficient(xi);
						// Compute the flux
						newFlux += ((double) size * coef * conc * hxLeft) / hy;
					}
				}
			}

			// Gather newFlux values at this position
			double newTotalFlux = 0.0;
			MPI_Allreduce(&newFlux, &newTotalFlux, 1, MPI_DOUBLE, MPI_SUM,
					PETSC_COMM_WORLD);

			// Update the previous flux
			previousIFlux2D[yj] = newTotalFlux;

			// Compare nInterstitials to the threshold to know if we should move the surface

			// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
			double threshold = (62.8 - initialVConc) * hxLeft * hy;
			if (nInterstitial2D[yj] > threshold) {
				// The surface is moving
				fvalue[0] = 0.0;
			}

			// Moving the surface back
			else if (nInterstitial2D[yj] < -threshold / 10.0) {
				// The surface is moving
				fvalue[0] = 0.0;
			}
		}
	}

	// Now work on the bubble bursting
	if (solverHandler.burstBubbles()) {
		// Compute the prefactor for the probability (arbitrary)
		double prefactor = heliumFluxAmplitude * dt * 0.1;

		// The depth parameter to know where the bursting should happen
		double depthParam = solverHandler.getTauBursting(); // nm

		// For now we are not bursting
		bool burst = false;

		// Loop on the full grid
		for (yj = 0; yj < My; yj++) {
			// Get the surface position
			int surfacePos = solverHandler.getSurfacePosition(yj);
			for (xi = 0; xi < Mx; xi++) {
				// Skip everything before the surface
				if (xi < surfacePos)
					continue;

				// If this is the locally owned part of the grid
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {

					// Get the pointer to the beginning of the solution data for this grid point
					gridPointSolution = solutionArray[yj][xi];
					// Update the concentration in the network
					network.updateConcentrationsFromArray(gridPointSolution);

					// Get the distance from the surface
					double distance = grid[xi + 1] - grid[surfacePos + 1];

					// Compute the helium density at this grid point
					double heDensity = network.getTotalAtomConcentration();

					// Compute the radius of the bubble from the number of helium
					double nV = heDensity * (grid[xi + 1] - grid[xi]) / 4.0;
					//				double nV = pow(heDensity / 5.0, 1.163) * (grid[xi + 1] - grid[xi]);
					double radius =
							(sqrt(3.0) / 4.0)
									* xolotlCore::tungstenLatticeConstant
									+ pow(
											(3.0
													* pow(
															xolotlCore::tungstenLatticeConstant,
															3.0) * nV)
													/ (8.0 * xolotlCore::pi),
											(1.0 / 3.0))
									- pow(
											(3.0
													* pow(
															xolotlCore::tungstenLatticeConstant,
															3.0))
													/ (8.0 * xolotlCore::pi),
											(1.0 / 3.0));

					// If the radius is larger than the distance to the surface, burst
					if (radius > distance) {
						burst = true;
						depthPositions2D.push_back(std::make_pair(yj, xi));
						// Exit the loop
						continue;
					}
					// Add randomness
					double prob = prefactor
							* (1.0 - (distance - radius) / distance)
							* min(1.0,
									exp(
											-(distance - depthParam)
													/ (depthParam * 2.0)));
					double test = solverHandler.getRNG().GetRandomDouble();

					if (prob > test) {
						burst = true;
						depthPositions2D.push_back(std::make_pair(yj, xi));
					}
				}
			}
		}

		// If at least one grid point is bursting
		if (burst) {
			// The event is happening
			fvalue[1] = 0.0;
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "postEventFunction2D")
/**
 * This is a method that moves the surface or burst bubbles
 */
PetscErrorCode postEventFunction2D(TS ts, PetscInt nevents,
		PetscInt eventList[], PetscReal time, Vec solution, PetscBool, void*) {

	// Initial declaration
	PetscErrorCode ierr;
	double ***solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi, Mx, ys, ym, yj, My;

	PetscFunctionBeginUser;

	// Check if the surface has moved
	if (nevents == 0)
		PetscFunctionReturn(0);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	CHKERRQ(ierr);

	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto& network = solverHandler.getNetwork();
	int dof = network.getDOF();

	// Get the physical grid
	auto grid = solverHandler.getXGrid();
	// Get the step size in Y
	double hy = solverHandler.getStepSizeY();

	// Take care of bursting

	// Loop on each bursting depth
	for (int i = 0; i < depthPositions2D.size(); i++) {
		// Get the coordinates of the point
		int xi = depthPositions2D[i].second, yj = depthPositions2D[i].first;
		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray[yj][xi];
		// Update the concentration in the network
		network.updateConcentrationsFromArray(gridPointSolution);

		// Get the surface position
		int surfacePos = solverHandler.getSurfacePosition(yj);
		// Get the distance from the surface
		double distance = grid[xi + 1] - grid[surfacePos + 1];

		std::cout << "bursting at: " << yj * hy << " " << distance << std::endl;

		// Pinhole case
		// Consider each He to reset their concentration at this grid point
		for (auto const& heMapItem : network.getAll(ReactantType::He)) {
			auto const& cluster = *(heMapItem.second);

			int id = cluster.getId() - 1;
			gridPointSolution[id] = 0.0;
		}
		// Consider each D to reset their concentration at this grid point
		for (auto const& dMapItem : network.getAll(ReactantType::D)) {
			auto const& cluster = *(dMapItem.second);

			int id = cluster.getId() - 1;
			gridPointSolution[id] = 0.0;
		}
		// Consider each T to reset their concentration at this grid point
		for (auto const& tMapItem : network.getAll(ReactantType::T)) {
			auto const& cluster = *(tMapItem.second);

			int id = cluster.getId() - 1;
			gridPointSolution[id] = 0.0;
		}

		// Consider each HeV cluster to transfer their concentration to the V cluster of the
		// same size at this grid point
		for (auto const& heVMapItem : network.getAll(ReactantType::PSIMixed)) {
			auto const& cluster = *(heVMapItem.second);

			// Get the V cluster of the same size
			auto const & comp = cluster.getComposition();
			auto vCluster = network.get(Species::V,
					comp[toCompIdx(Species::V)]);
			int vId = vCluster->getId() - 1;
			int id = cluster.getId() - 1;
			gridPointSolution[vId] += gridPointSolution[id];
			gridPointSolution[id] = 0.0;
		}

		// Loop on the super clusters to transfer their concentration to the V cluster of the
		// same size at this grid point
		for (auto const& superMapItem : network.getAll(ReactantType::PSISuper)) {
			auto const& cluster =
					static_cast<PSISuperCluster&>(*(superMapItem.second));

			// Loop on the V boundaries
			for (auto const& j : cluster.getBounds(3)) {
				// Get the total concentration at this v
				double conc = cluster.getIntegratedVConcentration(j);
				// Get the corresponding V cluster and its Id
				auto vCluster = network.get(Species::V, j);
				int vId = vCluster->getId() - 1;
				// Add the concentration
				gridPointSolution[vId] += conc;
			}

			// Reset the super cluster concentration
			int id = cluster.getId() - 1;
			gridPointSolution[id] = 0.0;
			id = cluster.getMomentId(0) - 1;
			gridPointSolution[id] = 0.0;
			id = cluster.getMomentId(1) - 1;
			gridPointSolution[id] = 0.0;
			id = cluster.getMomentId(2) - 1;
			gridPointSolution[id] = 0.0;
			id = cluster.getMomentId(3) - 1;
			gridPointSolution[id] = 0.0;
		}
	}

	// Now takes care of moving surface
	bool moving = false;
	for (int i = 0; i < nevents; i++) {
		if (eventList[i] == 0)
			moving = true;
	}

	// Skip if nothing is moving
	if (!moving) {
		// Restore the solutionArray
		ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
		CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}

	// Get the initial vacancy concentration
	double initialVConc = solverHandler.getInitialVConc();

	// Loop on the possible yj
	for (yj = 0; yj < My; yj++) {
		// Get the position of the surface at yj
		int surfacePos = solverHandler.getSurfacePosition(yj);
		xi = surfacePos + 1;

		// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
		double threshold = (62.8 - initialVConc) * (grid[xi + 1] - grid[xi])
				* hy;

		// Move the surface up
		if (nInterstitial2D[yj] > threshold) {
			int nGridPoints = 0;
			// Move the surface up until it is smaller than the next threshold
			while (nInterstitial2D[yj] > threshold) {
				// Move the surface higher
				surfacePos--;
				xi = surfacePos + 1;
				nGridPoints++;
				// Update the number of interstitials
				nInterstitial2D[yj] -= threshold;
				// Update the thresold
				double threshold = (62.8 - initialVConc)
						* (grid[xi + 1] - grid[xi]) * hy;
			}

			// Throw an exception if the position is negative
			if (surfacePos < 0) {
				throw std::string(
						"\nxolotlSolver::Monitor2D: The surface is trying to go outside of the grid!!");
			}

			// Printing information about the extension of the material
			if (procId == 0) {
				std::cout << "Adding " << nGridPoints
						<< " points to the grid on " << yj * hy << " at time: "
						<< time << " s." << std::endl;
			}

			// Set it in the solver
			solverHandler.setSurfacePosition(surfacePos, yj);

			// Initialize the vacancy concentration and the temperature on the new grid points
			// Get the single vacancy ID
			auto singleVacancyCluster = network.get(Species::V, 1);
			int vacancyIndex = -1;
			if (singleVacancyCluster)
				vacancyIndex = singleVacancyCluster->getId() - 1;
			// Get the surface temperature
			double temp = 0.0;
			if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
				temp = solutionArray[yj][xi][dof - 1];
			}
			double surfTemp = 0.0;
			MPI_Allreduce(&temp, &surfTemp, 1, MPI_DOUBLE, MPI_SUM,
					PETSC_COMM_WORLD);
			// Loop on the new grid points
			while (nGridPoints >= 0) {
				// Position of the newly created grid point
				xi = surfacePos + nGridPoints;

				// If xi is on this process
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
					// Get the concentrations
					gridPointSolution = solutionArray[yj][xi];

					// Set the new surface temperature
					gridPointSolution[dof - 1] = surfTemp;

					if (vacancyIndex > 0 && nGridPoints > 0) {
						// Initialize the vacancy concentration
						gridPointSolution[vacancyIndex] = initialVConc;
					}
				}

				// Decrease the number of grid points
				--nGridPoints;
			}
		}

		// Moving the surface back
		else if (nInterstitial2D[yj] < -threshold / 10.0) {
			// Move it back as long as the number of interstitials in negative
			while (nInterstitial2D[yj] < 0.0) {
				// Compute the threshold to a deeper grid point
				threshold = (62.8 - initialVConc)
						* (grid[xi + 2] - grid[xi + 1]) * hy;
				// Set all the concentrations to 0.0 at xi = surfacePos + 1
				// if xi is on this process
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
					// Get the concentrations at xi = surfacePos + 1
					gridPointSolution = solutionArray[yj][xi];
					// Loop on DOF
					for (int i = 0; i < dof - 1; i++) {
						gridPointSolution[i] = 0.0;
					}
				}

				// Move the surface deeper
				surfacePos++;
				xi = surfacePos + 1;
				// Update the number of interstitials
				nInterstitial2D[yj] += threshold;
			}

			// Printing information about the extension of the material
			if (procId == 0) {
				std::cout << "Removing grid points to the grid on " << yj * hy
						<< " at time: " << time << " s." << std::endl;
			}

			// Set it in the solver
			solverHandler.setSurfacePosition(surfacePos, yj);
		}
	}

	// Get the modified trap-mutation handler to reinitialize it
	auto mutationHandler = solverHandler.getMutationHandler();
	auto advecHandlers = solverHandler.getAdvectionHandlers();

	// Get the vector of positions of the surface
	std::vector<int> surfaceIndices;
	for (PetscInt i = 0; i < My; i++) {
		surfaceIndices.push_back(solverHandler.getSurfacePosition(i));
	}

	mutationHandler->initializeIndex2D(surfaceIndices, network, advecHandlers,
			grid, My, hy);

	// Write the surface positions
	if (procId == 0) {
		std::ofstream outputFile;
		outputFile.open("surface.txt", ios::app);
		outputFile << time << " ";

		// Loop on the possible yj
		for (yj = 0; yj < My; yj++) {
			// Get the position of the surface at yj
			int surfacePos = solverHandler.getSurfacePosition(yj);
			outputFile << grid[surfacePos + 1] - grid[1] << " ";
		}
		outputFile << std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/**
 * This operation sets up different monitors
 *  depending on the options.
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
PetscErrorCode setupPetsc2DMonitor(TS ts) {
	PetscErrorCode ierr;

	// Get the process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the xolotlViz handler registry
	auto vizHandlerRegistry = xolotlFactory::getVizHandlerRegistry();

	// Flags to launch the monitors or not
	PetscBool flagCheck, flagPerf, flagRetention, flagStatus, flag2DPlot,
			flagTRIDYN;

	// Check the option -check_collapse
	ierr = PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCheck);
	checkPetscError(ierr,
			"setupPetsc2DMonitor: PetscOptionsHasName (-check_collapse) failed.");

	// Check the option -plot_perf
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_perf", &flagPerf);
	checkPetscError(ierr,
			"setupPetsc2DMonitor: PetscOptionsHasName (-plot_perf) failed.");

	// Check the option -plot_2d
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_2d", &flag2DPlot);
	checkPetscError(ierr,
			"setupPetsc2DMonitor: PetscOptionsHasName (-plot_2d) failed.");

	// Check the option -helium_retention
	ierr = PetscOptionsHasName(NULL, NULL, "-helium_retention", &flagRetention);
	checkPetscError(ierr,
			"setupPetsc2DMonitor: PetscOptionsHasName (-helium_retention) failed.");

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, NULL, "-start_stop", &flagStatus);
	checkPetscError(ierr,
			"setupPetsc2DMonitor: PetscOptionsHasName (-start_stop) failed.");

	// Check the option -tridyn
	ierr = PetscOptionsHasName(NULL, NULL, "-tridyn", &flagTRIDYN);
	checkPetscError(ierr,
			"setupPetsc2DMonitor: PetscOptionsHasName (-tridyn) failed.");

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto& network = solverHandler.getNetwork();

	// Determine if we have an existing restart file,
	// and if so, it it has had timesteps written to it.
	std::unique_ptr<xolotlCore::XFile> networkFile;
	std::unique_ptr<xolotlCore::XFile::TimestepGroup> lastTsGroup;
	std::string networkName = solverHandler.getNetworkName();
	bool hasConcentrations = false;
	if (not networkName.empty()) {
		networkFile.reset(new xolotlCore::XFile(networkName));
		auto concGroup = networkFile->getGroup<
				xolotlCore::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
		if (hasConcentrations) {
			lastTsGroup = concGroup->getLastTimestepGroup();
		}
	}

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);
	checkPetscError(ierr, "setupPetsc2DMonitor: TSGetDM failed.");

	// Get the total size of the grid
	PetscInt Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);
	checkPetscError(ierr, "setupPetsc2DMonitor: DMDAGetInfo failed.");

	// Set the post step processing to stop the solver if the time step collapses
	if (flagCheck) {
		// Find the threshold
		PetscBool flag;
		ierr = PetscOptionsGetReal(NULL, NULL, "-check_collapse",
				&timeStepThreshold, &flag);
		checkPetscError(ierr,
				"setupPetsc2DMonitor: PetscOptionsGetInt (-check_collapse) failed.");
		if (!flag)
			timeStepThreshold = 1.0e-16;

		// Set the post step process that tells the solver when to stop if the time step collapse
		ierr = TSSetPostStep(ts, checkTimeStep);
		checkPetscError(ierr,
				"setupPetsc2DMonitor: TSSetPostStep (checkTimeStep) failed.");
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr = PetscOptionsGetReal(NULL, NULL, "-start_stop", &hdf5Stride2D,
				&flag);
		checkPetscError(ierr,
				"setupPetsc2DMonitor: PetscOptionsGetInt (-start_stop) failed.");
		if (!flag)
			hdf5Stride2D = 1.0;

		if (hasConcentrations) {
			// Get the previous time from the HDF5 file
			previousTime = lastTsGroup->readPreviousTime();
			hdf5Previous2D = (int) (previousTime / hdf5Stride2D);
		}

		// Don't do anything if both files have the same name
		if (hdf5OutputName2D != solverHandler.getNetworkName()) {

			// Get the solver handler
			auto& solverHandler = PetscSolver::getSolverHandler();

			// Get the physical grid in the x direction
			auto grid = solverHandler.getXGrid();

			// Setup step size variables
			double hy = solverHandler.getStepSizeY();

			// Get the compostion list and save it
			auto compList = network.getCompositionList();

			// Create and initialize a checkpoint file.
			// We do this in its own scope so that the file
			// is closed when the file object goes out of scope.
			// We want it to close before we (potentially) copy
			// the network from another file using a single-process
			// MPI communicator.
			{
				xolotlCore::XFile checkpointFile(hdf5OutputName2D, grid,
						compList, PETSC_COMM_WORLD, My, hy);
			}

			// Copy the network group from the given file (if it has one).
			// We open the files using a single-process MPI communicator
			// because it is faster for a single process to do the
			// copy with HDF5's H5Ocopy implementation than it is
			// when all processes call the copy function.
			// The checkpoint file must be closed before doing this.
			writeNetwork(PETSC_COMM_WORLD, solverHandler.getNetworkName(),
					hdf5OutputName2D, network);
		}

		// startStop2D will be called at each timestep
		ierr = TSMonitorSet(ts, startStop2D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc2DMonitor: TSMonitorSet (startStop2D) failed.");
	}

	// If the user wants the surface to be able to move or bursting
	if (solverHandler.moveSurface() || solverHandler.burstBubbles()) {
		// Surface
		if (solverHandler.moveSurface()) {
			// Initialize nInterstitial2D and previousIFlux2D before monitoring the
			// interstitial flux
			for (PetscInt j = 0; j < My; j++) {
				nInterstitial2D.push_back(0.0);
				previousIFlux2D.push_back(0.0);
			}

			// Get the interstitial information at the surface if concentrations were stored
			if (hasConcentrations) {

				assert(lastTsGroup);

				// Get the interstitial quantity from the HDF5 file
				nInterstitial2D = lastTsGroup->readData2D("nInterstitial");
				// Get the previous I flux from the HDF5 file
				previousIFlux2D = lastTsGroup->readData2D("previousIFlux");
				// Get the previous time from the HDF5 file
				previousTime = lastTsGroup->readPreviousTime();
			}

			// Get the sputtering yield
			sputteringYield2D = solverHandler.getSputteringYield();

			// Clear the file where the surface will be written
			std::ofstream outputFile;
			outputFile.open("surface.txt");
			outputFile.close();
		}

		// Bursting
		if (solverHandler.burstBubbles()) {
			// No need to seed the random number generator here.
			// The solver handler has already done it.
		}

		// Set directions and terminate flags for the surface event
		PetscInt direction[2];
		PetscBool terminate[2];
		direction[0] = 0, direction[1] = 0;
		terminate[0] = PETSC_FALSE, terminate[1] = PETSC_FALSE;
		// Set the TSEvent
		ierr = TSSetEventHandler(ts, 2, direction, terminate, eventFunction2D,
				postEventFunction2D, NULL);
		checkPetscError(ierr,
				"setupPetsc2DMonitor: TSSetEventHandler (eventFunction2D) failed.");
	}

	// Set the monitor to save performance plots (has to be in parallel)
	if (flagPerf) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a ScatterPlot
			perfPlot = vizHandlerRegistry->getPlot("perfPlot",
					xolotlViz::PlotType::SCATTER);

			// Create and set the label provider
			auto labelProvider = std::make_shared<xolotlViz::LabelProvider>(
					"labelProvider");
			labelProvider->axis1Label = "Process ID";
			labelProvider->axis2Label = "Solver Time";

			// Give it to the plot
			perfPlot->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider = std::make_shared<xolotlViz::CvsXDataProvider>(
					"dataProvider");

			// Give it to the plot
			perfPlot->setDataProvider(dataProvider);
		}

		// monitorPerf will be called at each timestep
		ierr = TSMonitorSet(ts, monitorPerf, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc2DMonitor: TSMonitorSet (monitorPerf) failed.");
	}

	// Set the monitor to compute the helium fluence for the retention calculation
	if (flagRetention) {
		// Check if we have a free surface at the bottom
		if (solverHandler.getRightOffset() == 1) {
			// Initialize n2D and previousFlux2D before monitoring the fluxes
			for (PetscInt j = 0; j < My; j++) {
				nHelium2D.push_back(0.0);
				previousHeFlux2D.push_back(0.0);
				nDeuterium2D.push_back(0.0);
				previousDFlux2D.push_back(0.0);
				nTritium2D.push_back(0.0);
				previousTFlux2D.push_back(0.0);
			}
		}

		// Get the previous time if concentrations were stored and initialize the fluence
		if (hasConcentrations) {

			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double time = lastTsGroup->readPreviousTime();
			// Initialize the fluence
			auto fluxHandler = solverHandler.getFluxHandler();
			// The length of the time step
			double dt = time;
			// Increment the fluence with the value at this current timestep
			fluxHandler->incrementFluence(dt);
			// Get the previous time from the HDF5 file
			// TODO is this the same as 'time' above?
			previousTime = lastTsGroup->readPreviousTime();

			// If the bottom is a free surface
			if (solverHandler.getRightOffset() == 1) {
				// Read about the impurity fluxes in the bulk
				nHelium2D = lastTsGroup->readData2D("nHelium");
				previousHeFlux2D = lastTsGroup->readData2D("previousHeFlux");
				nDeuterium2D = lastTsGroup->readData2D("nDeuterium");
				previousDFlux2D = lastTsGroup->readData2D("previousDFlux");
				nTritium2D = lastTsGroup->readData2D("nTritium");
				previousTFlux2D = lastTsGroup->readData2D("previousTFlux");
			}
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeFluence, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc2DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeHeliumRetention2D will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumRetention2D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc2DMonitor: TSMonitorSet (computeHeliumRetention2D) failed.");

		// Uncomment to clear the file where the retention will be written
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt");
		outputFile.close();
	}

	// Set the monitor to save surface plots of clusters concentration
	if (flag2DPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a SurfacePlot
			surfacePlot2D = vizHandlerRegistry->getPlot("surfacePlot2D",
					xolotlViz::PlotType::SURFACE);

			// Create and set the label provider
			auto labelProvider = std::make_shared<xolotlViz::LabelProvider>(
					"labelProvider");
			labelProvider->axis1Label = "Depth (nm)";
			labelProvider->axis2Label = "Y (nm)";
			labelProvider->axis3Label = "Concentration";

			// Give it to the plot
			surfacePlot2D->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider = std::make_shared<xolotlViz::CvsXYDataProvider>(
					"dataProvider");

			// Give it to the plot
			surfacePlot2D->setDataProvider(dataProvider);
		}

		// monitorSurface2D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSurface2D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc2DMonitor: TSMonitorSet (monitorSurface2D) failed.");
	}

	// Set the monitor to output data for TRIDYN
	if (flagTRIDYN) {
		// computeTRIDYN2D will be called at each timestep
		ierr = TSMonitorSet(ts, computeTRIDYN2D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc2DMonitor: TSMonitorSet (computeTRIDYN2D) failed.");
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(ts, monitorTime, NULL, NULL);
	checkPetscError(ierr,
			"setupPetsc2DMonitor: TSMonitorSet (monitorTime) failed.");

	PetscFunctionReturn(0);
}

}

/* end namespace xolotlSolver */
