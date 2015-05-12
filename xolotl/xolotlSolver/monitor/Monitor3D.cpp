// Includes
#include "PetscSolver.h"
#include <xolotlPerf.h>
#include <VizHandlerRegistryFactory.h>
#include <PlotType.h>
#include <CvsXDataProvider.h>
#include <CvsXYDataProvider.h>
#include <LabelProvider.h>
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <memory>
#include <HDF5Utils.h>

namespace xolotlSolver {

// Declaration of the functions defined in Monitor.cpp
extern PetscErrorCode monitorTime(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode computeHeliumFluence(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode monitorPerf(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);

// Declaration of the variables defined in Monitor.cpp
extern std::shared_ptr<xolotlViz::IPlot> perfPlot;
extern double previousTime;

//! How often HDF5 file is written
PetscInt hdf5Stride3D = 0;
//! HDF5 output file name
std::string hdf5OutputName3D = "xolotlStop.h5";
// Declare the vector that will store the Id of the helium clusters
std::vector<int> heIndices3D;
// Declare the vector that will store the weight of the helium clusters
// (their He composition)
std::vector<int> heWeights3D;
//! The pointer to the 2D plot used in MonitorSurfaceXY3D.
std::shared_ptr<xolotlViz::IPlot> surfacePlotXY3D;
//! The pointer to the 2D plot used in MonitorSurfaceXZ3D.
std::shared_ptr<xolotlViz::IPlot> surfacePlotXZ3D;

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "startStop3D")
/**
 * This is a monitoring method that will save an hdf5 file at each time step.
 * HDF5 is handling the parallel part, so no call to MPI here.
 */
PetscErrorCode startStop3D(TS ts, PetscInt timestep, PetscReal time, Vec solution,
		void *ictx) {
	PetscErrorCode ierr;
	const double ****solutionArray, *gridPointSolution;
	int xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % hdf5Stride3D != 0)
		PetscFunctionReturn(0);

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);CHKERRQ(ierr);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto network = solverHandler->getNetwork();

	// Network size
	const int networkSize = network->size();

	// Open the already created HDF5 file
	xolotlCore::HDF5Utils::openFile(hdf5OutputName3D);

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);CHKERRQ(ierr);

	// Add a concentration sub group
	xolotlCore::HDF5Utils::addConcentrationSubGroup(timestep, networkSize, time,
			currentTimeStep);

	// Loop on the full grid
	for (int k = 0; k < Mz; k++) {
		for (int j = 0; j < My; j++) {
			for (int i = 0; i < Mx; i++) {
				// Size of the concentration that will be stored
				int concSize = -1;
				// Vector for the concentrations
				std::vector<std::vector<double> > concVector;

				// If it is the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym
						&& k >= zs && k < zs + zm) {
					// Get the pointer to the beginning of the solution data
					// for this grid point
					gridPointSolution = solutionArray[k][j][i];

					// Loop on the concentrations
					concVector.clear();
					for (int l = 0; l < networkSize; l++) {
						if (gridPointSolution[l] > 1.0e-16) {
							// Create the concentration vector for this cluster
							std::vector<double> conc;
							conc.push_back((double) l);
							conc.push_back(gridPointSolution[l]);

							// Add it to the main vector
							concVector.push_back(conc);
						}
					}

					// Send the size of the vector to the other processes
					concSize = concVector.size();
					// Loop on all the processes
					for (int l = 0; l < worldSize; l++) {
						// Skip its own
						if (l == procId)
							continue;

						// Send the size
						MPI_Send(&concSize, 1, MPI_INT, l, 0, MPI_COMM_WORLD);
					}
				}

				// Else: only receive the conc size
				else {
					MPI_Recv(&concSize, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
							MPI_STATUS_IGNORE);
				}

				// Skip the grid point if the size is 0
				if (concSize == 0)
					continue;

				// All processes must create the dataset
				xolotlCore::HDF5Utils::addConcentrationDataset(concSize, i, j, k);

				// If it is the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym
						&& k >= zs && k < zs + zm) {
					// Fill the dataset
					xolotlCore::HDF5Utils::fillConcentrations(concVector, i, j, k);
				}
			}
		}
	}

	// Finalize the HDF5 file
	xolotlCore::HDF5Utils::closeFile();

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeHeliumRetention3D")
/**
 * This is a monitoring method that will compute the total helium fluence
 */
PetscErrorCode computeHeliumRetention3D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	PetscErrorCode ierr;
	int xs, xm, ys, ym, zs, zm;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = solverHandler->getFluxHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	// Setup step size variables
	double hx = solverHandler->getStepSizeX();
	double hy = solverHandler->getStepSizeY();
	double hz = solverHandler->getStepSizeZ();

	// Get the array of concentration
	const double ****solutionArray, *gridPointSolution;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	// Store the concentration over the grid
	double heConcentration = 0;

	// Loop on the grid
	for (int k = zs; k < zs + zm; k++) {
		for (int j = ys; j < ys + ym; j++) {
			for (int i = xs; i < xs + xm; i++) {
				// Get the pointer to the beginning of the solution data for
				// this grid point
				gridPointSolution = solutionArray[k][j][i];

				// Loop on all the indices
				for (int l = 0; l < heIndices3D.size(); l++) {
					// Add the current concentration times the number of
					// helium in the cluster (from the weight vector)
					heConcentration += gridPointSolution[heIndices3D[l]]
					                                     * heWeights3D[l] * hx;
				}
			}
		}
	}

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);
	// Get the current process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Sum all the concentrations through MPI reduce
	double totalHeConcentration = 0.0;
	MPI_Reduce(&heConcentration, &totalHeConcentration, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// Master process
	if (procId == 0) {

		// Get the total size of the grid rescale the concentrations
		int Mx, My, Mz;
		ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz,
				PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
				PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
				PETSC_IGNORE);CHKERRQ(ierr);

		// Compute the total surface irradiated by the helium flux
		double surface = (double) (My * Mz) * hy * hz;

		// Rescale the concentration
		totalHeConcentration = totalHeConcentration / surface;

		// Get the fluence
		double heliumFluence = fluxHandler->getHeFluence();

		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		std::cout << "Helium retention = "
				<< 100.0 * (totalHeConcentration / heliumFluence) << " %"
				<< std::endl;
		std::cout << "Helium mean concentration = " << totalHeConcentration << std::endl;
		std::cout << "Helium fluence = " << heliumFluence << "\n" << std::endl;

//		// Uncomment to write the retention and the fluence in a file
//		std::ofstream outputFile;
//		outputFile.open("retentionOut.txt", ios::app);
//		outputFile << heliumFluence << " "
//				<< 100.0 * (totalHeConcentration / heliumFluence) << std::endl;
//		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorSurfaceXY3D")
/**
 * This is a monitoring method that will save 2D plots of the concentration of
 * a specific cluster at each grid point on the XY surface, integrating over Z.
 */
PetscErrorCode monitorSurfaceXY3D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	PetscErrorCode ierr;
	const double ****solutionArray, *gridPointSolution;
	int xs, xm, Mx, ys, ym, My, zs, zm, Mz;
	double x, y;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);CHKERRQ(ierr);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto network = solverHandler->getNetwork();

	// Setup step size variables
	double hx = solverHandler->getStepSizeX();
	double hy = solverHandler->getStepSizeY();
	double hz = solverHandler->getStepSizeZ();

	// Choice of the cluster to be plotted
	int iCluster = 0;

	// Create a Point vector to store the data to give to the data provider
	// for the visualization
	auto myPoints = std::make_shared<std::vector<xolotlViz::Point> >();
	// Create a point here so that it is not created and deleted in the loop
	xolotlViz::Point thePoint;

	// Loop on the full grid, Y and X first because they are the axis of the plot
	for (int j = 0; j < My; j++) {
		// Compute y
		y = (double) j * hy;

		for (int i = 0; i < Mx; i++) {
			// Compute x
			x = (double) i * hx;

			// Initialize the value of the concentration to integrate over Z
			double conc = 0.0;

			for (int k = 0; k < Mz; k++) {
				// If it is the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym
						&& k >= zs && k < zs + zm) {
					// Get the pointer to the beginning of the solution data for this grid point
					gridPointSolution = solutionArray[k][j][i];

					// Integrate over Z
					conc += gridPointSolution[iCluster];
				}
			} // End of the loop on Z

			// Sum all the concentration on Z
			double totalConc = 0.0;
			MPI_Reduce(&conc, &totalConc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

			// If it is procId == 0 just store the integrated value in the myPoints vector
			if (procId == 0) {
				thePoint.value = totalConc;
				thePoint.t = time;
				thePoint.x = x;
				thePoint.y = y;
				myPoints->push_back(thePoint);
			}
		}
	}

	// Plot everything from procId == 0
	if (procId == 0) {
		// Get the data provider and give it the points
		surfacePlotXY3D->getDataProvider()->setPoints(myPoints);

		// Get the iCluster cluster to have access to its name
		auto reactants = network->getAll();
		auto cluster = (PSICluster *) reactants->at(iCluster);

		// Change the title of the plot and the name of the data
		std::stringstream title;
		title << cluster->getName();
		surfacePlotXY3D->getDataProvider()->setDataName(title.str());
		title << " concentration";
		surfacePlotXY3D->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		surfacePlotXY3D->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		surfacePlotXY3D->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << cluster->getName() << "_surfaceXY_TS" << timestep << ".pnm";
		surfacePlotXY3D->write(fileName.str());
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorSurfaceXZ3D")
/**
 * This is a monitoring method that will save 2D plots of the concentration of
 * a specific cluster at each grid point on the XZ surface, integrating over Y.
 */
PetscErrorCode monitorSurfaceXZ3D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	PetscErrorCode ierr;
	const double ****solutionArray, *gridPointSolution;
	int xs, xm, Mx, ys, ym, My, zs, zm, Mz;
	double x, z;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);CHKERRQ(ierr);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto network = solverHandler->getNetwork();

	// Setup step size variables
	double hx = solverHandler->getStepSizeX();
	double hy = solverHandler->getStepSizeY();
	double hz = solverHandler->getStepSizeZ();

	// Choice of the cluster to be plotted
	int iCluster = 0;

	// Create a Point vector to store the data to give to the data provider
	// for the visualization
	auto myPoints = std::make_shared<std::vector<xolotlViz::Point> >();
	// Create a point here so that it is not created and deleted in the loop
	xolotlViz::Point thePoint;

	// Loop on the full grid, Z and X first because they are the axis of the plot
	for (int k = 0; k < Mz; k++) {
		// Compute z
		z = k * hz;

		for (int i = 0; i < Mx; i++) {
			// Compute x
			x = i * hx;

			// Initialize the value of the concentration to integrate over Y
			double conc = 0.0;

			for (int j = 0; j < My; j++) {
				// If it is the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym
						&& k >= zs && k < zs + zm) {
					// Get the pointer to the beginning of the solution data for this grid point
					gridPointSolution = solutionArray[k][j][i];

					// Integrate over Y
					conc += gridPointSolution[iCluster];
				}
			} // End of the loop on Y

			// Sum all the concentration on Y
			double totalConc = 0.0;
			MPI_Reduce(&conc, &totalConc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

			// If it is procId == 0 just store the integrated value in the myPoints vector
			if (procId == 0) {
				thePoint.value = totalConc;
				thePoint.t = time;
				thePoint.x = x;
				thePoint.y = z;
				myPoints->push_back(thePoint);
			}
		}
	}

	// Plot everything from procId == 0
	if (procId == 0) {
		// Get the data provider and give it the points
		surfacePlotXZ3D->getDataProvider()->setPoints(myPoints);

		// Get the iCluster cluster to have access to its name
		auto reactants = network->getAll();
		auto cluster = (PSICluster *) reactants->at(iCluster);

		// Change the title of the plot and the name of the data
		std::stringstream title;
		title << cluster->getName();
		surfacePlotXZ3D->getDataProvider()->setDataName(title.str());
		title << " concentration";
		surfacePlotXZ3D->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		surfacePlotXZ3D->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		surfacePlotXZ3D->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << cluster->getName() << "_surfaceXZ_TS" << timestep << ".pnm";
		surfacePlotXZ3D->write(fileName.str());
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/**
 * This operation sets up all the needed monitors.
 *
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
PetscErrorCode setupPetsc3DMonitor(TS ts) {
	PetscErrorCode ierr;

	// Get the process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the xolotlViz handler registry
	auto vizHandlerRegistry = xolotlFactory::getVizHandlerRegistry();

	// Flags to launch the monitors or not
	PetscBool flagPerf, flagRetention, flagStatus, flag2DXYPlot, flag2DXZPlot;

	// Check the option -plot_perf
	ierr = PetscOptionsHasName(NULL, "-plot_perf", &flagPerf);CHKERRQ(ierr);

	// Check the option -plot_2d_xy
	ierr = PetscOptionsHasName(NULL, "-plot_2d_xy", &flag2DXYPlot);CHKERRQ(ierr);

	// Check the option -plot_2d_xz
	ierr = PetscOptionsHasName(NULL, "-plot_2d_xz", &flag2DXZPlot);CHKERRQ(ierr);

	// Check the option -helium_retention
	ierr = PetscOptionsHasName(NULL, "-helium_retention", &flagRetention);CHKERRQ(ierr);

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, "-start_stop", &flagStatus);CHKERRQ(ierr);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto network = solverHandler->getNetwork();
	const int networkSize = network->size();

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
		ierr = TSMonitorSet(ts, monitorPerf, NULL, NULL);CHKERRQ(ierr);
	}

	// Set the monitor to compute the helium fluence for the retention calculation
	if (flagRetention) {
		// Get all the helium clusters
		auto heClusters = network->getAll(heType);

		// Get all the helium-vacancy clusters
		auto heVClusters = network->getAll(heVType);

		// Loop on the helium clusters
		for (int i = 0; i < heClusters.size(); i++) {
			auto cluster = (PSICluster *) heClusters[i];
			int id = cluster->getId() - 1;
			// Add the Id to the vector
			heIndices3D.push_back(id);
			// Add the number of heliums of this cluster to the weight
			heWeights3D.push_back(cluster->getSize());
		}

		// Loop on the helium-vacancy clusters
		for (int i = 0; i < heVClusters.size(); i++) {
			auto cluster = (PSICluster *) heVClusters[i];
			int id = cluster->getId() - 1;
			// Add the Id to the vector
			heIndices3D.push_back(id);
			// Add the number of heliums of this cluster to the weight
			auto comp = cluster->getComposition();
			heWeights3D.push_back(comp[heType]);
		}

		if (heIndices3D.size() == 0) {
			throw std::string(
					"PetscSolver Exception: Cannot compute the retention because "
					"there is no helium or helium-vacancy cluster in the network.");
		}

		// computeHeliumFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumFluence, NULL, NULL);CHKERRQ(ierr);

		// computeHeliumRetention3D will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumRetention3D, NULL, NULL);CHKERRQ(ierr);

//		// Uncomment to clear the file where the retention will be written
//		std::ofstream outputFile;
//		outputFile.open("retentionOut.txt");
//		outputFile.close();
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr = PetscOptionsGetInt(NULL, "-start_stop", &hdf5Stride3D, &flag);CHKERRQ(ierr);
		if (!flag)
			hdf5Stride3D = 1;

		PetscInt Mx, My, Mz;
		PetscErrorCode ierr;

		// Get the da from ts
		DM da;
		ierr = TSGetDM(ts, &da);CHKERRQ(ierr);

		// Get the size of the total grid
		ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz,
				PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
				PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
				PETSC_IGNORE);CHKERRQ(ierr);

		// Initialize the HDF5 file for all the processes
		xolotlCore::HDF5Utils::initializeFile(hdf5OutputName3D, networkSize);

		// Get the solver handler
		auto solverHandler = PetscSolver::getSolverHandler();

		// Setup step size variables
		double hx = solverHandler->getStepSizeX();
		double hy = solverHandler->getStepSizeY();
		double hz = solverHandler->getStepSizeZ();

		// Save the header in the HDF5 file
		xolotlCore::HDF5Utils::fillHeader(Mx, hx, My, hy, Mz, hz);

		// Save the network in the HDF5 file
		xolotlCore::HDF5Utils::fillNetwork(network);

		// Finalize the HDF5 file
		xolotlCore::HDF5Utils::finalizeFile();

		// startStop3D will be called at each timestep
		ierr = TSMonitorSet(ts, startStop3D, NULL, NULL);CHKERRQ(ierr);
	}

	// Set the monitor to save surface plots of clusters concentration
	if (flag2DXYPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a SurfacePlot
			surfacePlotXY3D = vizHandlerRegistry->getPlot("surfacePlotXY3D",
					xolotlViz::PlotType::SURFACE);

			// Create and set the label provider
			auto labelProvider = std::make_shared<xolotlViz::LabelProvider>(
					"labelProvider");
			labelProvider->axis1Label = "Depth (nm)";
			labelProvider->axis2Label = "Y (nm)";
			labelProvider->axis3Label = "Concentration";

			// Give it to the plot
			surfacePlotXY3D->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider = std::make_shared<xolotlViz::CvsXYDataProvider>(
					"dataProvider");

			// Give it to the plot
			surfacePlotXY3D->setDataProvider(dataProvider);
		}

		// monitorSurfaceXY3D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSurfaceXY3D, NULL, NULL);CHKERRQ(ierr);
	}

	// Set the monitor to save surface plots of clusters concentration
	if (flag2DXZPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a SurfacePlot
			surfacePlotXZ3D = vizHandlerRegistry->getPlot("surfacePlotXZ3D",
					xolotlViz::PlotType::SURFACE);

			// Create and set the label provider
			auto labelProvider = std::make_shared<xolotlViz::LabelProvider>(
					"labelProvider");
			labelProvider->axis1Label = "Depth (nm)";
			labelProvider->axis2Label = "Z (nm)";
			labelProvider->axis3Label = "Concentration";

			// Give it to the plot
			surfacePlotXZ3D->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider = std::make_shared<xolotlViz::CvsXYDataProvider>(
					"dataProvider");

			// Give it to the plot
			surfacePlotXZ3D->setDataProvider(dataProvider);
		}

		// monitorSurfaceXZ3D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSurfaceXZ3D, NULL, NULL);CHKERRQ(ierr);
	}

	// Set the monitor to simply change the previous time to the new time
	if (flagRetention) {
		// monitorTime will be called at each timestep
		ierr = TSMonitorSet(ts, monitorTime, NULL, NULL);CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}

}

/* end namespace xolotlSolver */
