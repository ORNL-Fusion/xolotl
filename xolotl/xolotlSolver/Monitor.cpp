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

/* ----- Error Handling Code ----- */

/**
 * This operation checks a Petsc error code and converts it to a bool.
 * @param errorCode The Petsc error code.
 * @return True if everything is OK, false otherwise.
 */
static inline bool checkPetscError(PetscErrorCode errorCode) {
	CHKERRQ(errorCode);
}

//! The pointer to the plot that will be used to visualize the data.
std::shared_ptr<xolotlViz::IPlot> plot;

//! The pointer to the series plot that will be used to visualize the data.
std::shared_ptr<xolotlViz::IPlot> seriesPlot;

//! The pointer to the 2D plot that will be used to visualize the data.
std::shared_ptr<xolotlViz::IPlot> surfacePlot;

//! The pointer to the plot that will be used to visualize performance data.
std::shared_ptr<xolotlViz::IPlot> perfPlot;

//! The variable to store the time at the previous time step.
double previousTime = 0.0;

//! The variable to store the interstitial flux at the previous time step.
double previousFlux = 0.0;

//! How often HDF5 file is written
PetscInt stride = 0;

//! HDF5 output file name
std::string outputFileName = "xolotlStop.h5";

// Declare the vector that will store the Id of the helium clusters
std::vector<int> indices;

// Declare the vector that will store the weight of the helium clusters
// (their He composition)
std::vector<int> weight;

// Variable to indicate whether or not the fact that the concentration of the biggest
// cluster in the network is higher than 1.0e-16 should be printed.
// Becomes false once it is printed.
bool printMaxClusterConc = true;

/**
 * This is a monitoring method set the previous time to the time. This is needed here
 * because multiple monitors need the previous time value from the previous timestep.
 * This monitor is called last.
 */
PetscErrorCode monitorTime(TS ts, PetscInt timestep, PetscReal time, Vec solution,
		void *ictx) {
	PetscFunctionBeginUser;

	// Set the previous time to the current time for the next timestep
	previousTime = time;

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will save an hdf5 file at each time step.
 * HDF5 is handling the parallel part, so no call to MPI here.
 */
PetscErrorCode startStop(TS ts, PetscInt timestep, PetscReal time, Vec solution,
		void *ictx) {
	// Network size
	const int networkSize = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	double *solutionArray, *gridPointSolution;
	Vec localSolution;
	int xs, xm, Mx;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % stride != 0)
		PetscFunctionReturn(0);

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

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

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Setup step size variable
	double hx = PetscSolver::getStepSize();

	// Open the already created HDF5 file
	xolotlCore::HDF5Utils::openFile(outputFileName);

	// Get the physical dimension of the grid
	int dimension = (Mx - 1) * hx;

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	checkPetscError(ierr);

	// Add a concentration sub group
	xolotlCore::HDF5Utils::addConcentrationSubGroup(timestep, networkSize, time,
			currentTimeStep);

	// Loop on the full grid
	for (int xi = 0; xi < Mx; xi++) {
		// Size of the concentration that will be stored
		int concSize = -1;
		// Vector for the concentrations
		std::vector<std::vector<double> > concVector;

		// If it is the locally owned part of the grid
		if (xi >= xs && xi < xs + xm) {
			// Get the pointer to the beginning of the solution data for this grid point
			gridPointSolution = solutionArray + networkSize * xi;
			// Update the concentrations in the network to have physics results
			// (non negative)
			PetscSolver::getNetwork()->updateConcentrationsFromArray(
					gridPointSolution);
			// Get the concentrations from the network
			double concentrations[networkSize];
			double * concentration = &concentrations[0];
			PetscSolver::getNetwork()->fillConcentrationsArray(concentration);

			// Loop on the concentrations
			concVector.clear();
			for (int i = 0; i < networkSize; i++) {
				if (concentrations[i] > 1.0e-16) {
					// Create the concentration vector for this cluster
					std::vector<double> conc;
					conc.push_back((double) i);
					conc.push_back(concentrations[i]);

					// Add it to the main vector
					concVector.push_back(conc);
				}
			}

			// Send the size of the vector to the other processes
			concSize = concVector.size();
			// Loop on all the processes
			for (int i = 0; i < worldSize; i++) {
				// Skip its own
				if (i == procId)
					continue;

				// Send the size
				MPI_Send(&concSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
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
		xolotlCore::HDF5Utils::addConcentrationDataset(xi, concSize);

		// If it is the locally owned part of the grid
		if (xi >= xs && xi < xs + xm) {
			// Fill the dataset
			xolotlCore::HDF5Utils::fillConcentrations(concVector, xi);
		}
	}

	// Finalize the HDF5 file
	xolotlCore::HDF5Utils::closeFile();

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will compute the total helium fluence
 */
PetscErrorCode computeHeliumFluence(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	//
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = PetscSolver::getFluxHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// The length of the time step
	double dt = time - previousTime;

	// Increment the fluence with the value at this current timestep
	fluxHandler->incrementHeFluence(dt);

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will compute the total helium fluence
 */
PetscErrorCode computeHeliumRetention(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int size = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	int xs, xm;

	PetscFunctionBeginUser;

	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = PetscSolver::getFluxHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Setup step size variable
	double hx = PetscSolver::getStepSize();

	// Get the array of concentration
	PetscReal *solutionArray;
	ierr = DMDAVecGetArray(da, solution, &solutionArray);
	checkPetscError(ierr);

	// Store the concentration over the grid
	double heConcentration = 0;

	// Loop on the grid
	for (int xi = xs; xi < xs + xm; xi++) {
		// Get the pointer to the beginning of the solution data for this grid point
		PetscReal *gridPointSolution;
		gridPointSolution = solutionArray + size * xi;

		// Update the concentrations in the network to have physics results
		// (non negative)
		PetscSolver::getNetwork()->updateConcentrationsFromArray(
				gridPointSolution);

		// Get the concentrations from the network
		double concentrations[size];
		double * concentration = &concentrations[0];
		PetscSolver::getNetwork()->fillConcentrationsArray(concentration);

		// Loop on all the indices
		for (int i = 0; i < indices.size(); i++) {
			// Add the current concentration times the number of helium in the cluster
			// (from the weight vector)
			heConcentration += concentration[indices[i]] * weight[i] * hx;
		}
	}

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);
	// Get the current process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Master process
	if (procId == 0) {
		// Loop on all the other processes
		for (int i = 1; i < worldSize; i++) {
			double otherConcentration = 0.0;

			// Receive the value from the other processes
			MPI_Recv(&otherConcentration, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			// Add them to the master one
			heConcentration += otherConcentration;
		}

		// Get the fluence
		double heliumFluence = fluxHandler->getHeFluence();

		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		std::cout << "Helium retention = "
				<< 100.0 * (heConcentration / heliumFluence) << " %"
				<< std::endl;
		std::cout << "Helium concentration = " << heConcentration << std::endl;
		std::cout << "Helium fluence = " << heliumFluence << "\n" << std::endl;

//		// Uncomment to write the retention and the fluence in a file
//		std::ofstream outputFile;
//		outputFile.open("retentionOut.txt", ios::app);
//		outputFile << heliumFluence << " "
//				<< 100.0 * (heConcentration / heliumFluence) << std::endl;
//		outputFile.close();
	}

	else {
		// Send the value of the timer to the master process
		MPI_Send(&heConcentration, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will save 1D plots of one concentration
 */
PetscErrorCode monitorScatter(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int networkSize = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	double *solutionArray, *gridPointSolution, x;
	Vec localSolution;
	int xs, xm, xi;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

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

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Setup some step size variables
	double hx = PetscSolver::getStepSize();

	// Choice of the cluster to be plotted
	int iCluster = 6;

	if (procId == 0) {
		// Create a Point vector to store the data to give to the data provider
		// for the visualization
		auto myPoints = std::make_shared<std::vector<xolotlViz::Point> >();

		// Loop on the grid
		for (xi = xs; xi < xs + xm; xi++) {
			// Dump x
			x = xi * hx;
			// Get the pointer to the beginning of the solution data for this grid point
			gridPointSolution = solutionArray + networkSize * xi;
			// Update the concentrations in the network to have physics results
			// (non negative)
			PetscSolver::getNetwork()->updateConcentrationsFromArray(
					gridPointSolution);
			// Get the concentrations from the network
			double concentrations[networkSize];
			double * concentration = &concentrations[0];
			PetscSolver::getNetwork()->fillConcentrationsArray(concentration);

			// Create a Point with the concentration[iCluster] as the value
			// and add it to myPoints
			xolotlViz::Point aPoint;
			aPoint.value = concentration[iCluster];
			aPoint.t = time;
			aPoint.x = x;
			myPoints->push_back(aPoint);
		}

		// Loop on the other processes
		for (int i = 1; i < worldSize; i++) {
			// Get the size of the local grid of that process
			int localSize = 0;
			MPI_Recv(&localSize, 1, MPI_INT, i, 2, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			// Loop on their grid
			for (int k = 0; k < localSize; k++) {
				// Get the position
				MPI_Recv(&x, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

				// and the concentration
				double conc = 0.0;
				MPI_Recv(&conc, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

				// Create a Point with the concentration[iCluster] as the value
				// and add it to myPoints
				xolotlViz::Point aPoint;
				aPoint.value = conc;
				aPoint.t = time;
				aPoint.x = x;
				myPoints->push_back(aPoint);
			}
		}

		// Get the data provider and give it the points
		plot->getDataProvider()->setPoints(myPoints);

		// Get the iCluster cluster to have access to its name
		auto reactants = PetscSolver::getNetwork()->getAll();
		auto cluster = (PSICluster *) reactants->at(iCluster);

		// Change the title of the plot and the name of the data
		std::stringstream title;
		title << cluster->getName();
		plot->getDataProvider()->setDataName(title.str());
		title << " concentration";
		plot->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		plot->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		checkPetscError(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		plot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << cluster->getName() << "_scatter_TS" << timestep << ".pnm";
		plot->write(fileName.str());
	}

	else {
		// Send the value of the local grid size to the master process
		MPI_Send(&xm, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);

		// Loop on the grid
		for (xi = xs; xi < xs + xm; xi++) {
			// Dump x
			x = xi * hx;

			// Get the pointer to the beginning of the solution data for this grid point
			gridPointSolution = solutionArray + networkSize * xi;
			// Update the concentrations in the network to have physics results
			// (non negative)
			PetscSolver::getNetwork()->updateConcentrationsFromArray(
					gridPointSolution);
			// Get the concentrations from the network
			double concentrations[networkSize];
			double * concentration = &concentrations[0];
			PetscSolver::getNetwork()->fillConcentrationsArray(concentration);

			// Send the value of the local position to the master process
			MPI_Send(&x, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);

			// Send the value of the concentration to the master process
			MPI_Send(&concentration[iCluster], 1, MPI_DOUBLE, 0, 2,
					MPI_COMM_WORLD);
		}
	}

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will save 1D plots of many concentrations
 */
PetscErrorCode monitorSeries(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int networkSize = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	double *solutionArray, *gridPointSolution, x;
	Vec localSolution;
	int xs, xm, xi;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

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

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Setup some step size variables
	double hx = PetscSolver::getStepSize();

	// To plot a maximum of 18 clusters of the whole benchmark
	const int loopSize = std::min(18, networkSize);

	if (procId == 0) {
		// Create a Point vector to store the data to give to the data provider
		// for the visualization
		std::vector<std::vector<xolotlViz::Point> > myPoints(loopSize);

		// Loop on the grid
		for (xi = xs; xi < xs + xm; xi++) {
			// Dump x
			x = xi * hx;
			// Get the pointer to the beginning of the solution data for this grid point
			gridPointSolution = solutionArray + networkSize * xi;
			// Update the concentrations in the network to have physics results
			// (non negative)
			PetscSolver::getNetwork()->updateConcentrationsFromArray(
					gridPointSolution);
			// Get the concentrations from the network
			double concentrations[networkSize];
			double * concentration = &concentrations[0];
			PetscSolver::getNetwork()->fillConcentrationsArray(concentration);

			for (int i = 0; i < loopSize; i++) {
				// Create a Point with the concentration[iCluster] as the value
				// and add it to myPoints
				xolotlViz::Point aPoint;
				aPoint.value = concentration[i];
				aPoint.t = time;
				aPoint.x = x;
				myPoints[i].push_back(aPoint);
			}
		}

		// Loop on the other processes
		for (int i = 1; i < worldSize; i++) {
			// Get the size of the local grid of that process
			int localSize = 0;
			MPI_Recv(&localSize, 1, MPI_INT, i, 3, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			// Loop on their grid
			for (int k = 0; k < localSize; k++) {
				// Get the position
				MPI_Recv(&x, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

				for (int j = 0; j < loopSize; j++) {
					// and the concentrations
					double conc;
					MPI_Recv(&conc, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD,
							MPI_STATUS_IGNORE);

					// Create a Point with the concentration[iCluster] as the value
					// and add it to myPoints
					xolotlViz::Point aPoint;
					aPoint.value = conc; // He
					aPoint.t = time;
					aPoint.x = x;
					myPoints[j].push_back(aPoint);
				}
			}
		}

		// Get all the reactants to have access to their names
		auto reactants = PetscSolver::getNetwork()->getAll();

		for (int i = 0; i < loopSize; i++) {
			auto cluster = (PSICluster *) reactants->at(i);
			// Get the data provider and give it the points
			auto thePoints = std::make_shared<std::vector<xolotlViz::Point> >(
					myPoints[i]);
			seriesPlot->getDataProvider(i)->setPoints(thePoints);
			seriesPlot->getDataProvider(i)->setDataName(cluster->getName());
		}

		// Change the title of the plot
		std::stringstream title;
		title << "Concentrations";
		seriesPlot->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		seriesPlot->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		checkPetscError(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		seriesPlot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "log_series_TS" << timestep << ".pnm";
		seriesPlot->write(fileName.str());
	}

	else {
		// Send the value of the local grid size to the master process
		MPI_Send(&xm, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

		// Loop on the grid
		for (xi = xs; xi < xs + xm; xi++) {
			// Dump x
			x = xi * hx;

			// Get the pointer to the beginning of the solution data for this grid point
			gridPointSolution = solutionArray + networkSize * xi;
			// Update the concentrations in the network to have physics results
			// (non negative)
			PetscSolver::getNetwork()->updateConcentrationsFromArray(
					gridPointSolution);
			// Get the concentrations from the network
			double concentrations[networkSize];
			double * concentration = &concentrations[0];
			PetscSolver::getNetwork()->fillConcentrationsArray(concentration);

			// Send the value of the local position to the master process
			MPI_Send(&x, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

			for (int i = 0; i < loopSize; i++) {
				// Send the value of the concentrations to the master process
				MPI_Send(&concentration[i], 1, MPI_DOUBLE, 0, 3,
						MPI_COMM_WORLD);
			}
		}
	}

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will save 2D plots for each depths of
 * the concentration as a function of the cluster composition.
 */
PetscErrorCode monitorSurface(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int networkSize = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	double *solutionArray, *gridPointSolution, x;
	Vec localSolution;
	int xs, xm, xi;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

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

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Setup some step size variables
	double hx = PetscSolver::getStepSize();

	// Get the maximum size of HeV clusters
	auto psiNetwork = std::dynamic_pointer_cast<PSIClusterReactionNetwork>(
			PetscSolver::getNetwork());
	std::map<std::string, std::string> props = psiNetwork->getProperties();
	int maxHeVClusterSize = std::stoi(props["maxHeVClusterSize"]);
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);

	// Loop on the grid points
	for (xi = xs; xi < xs + xm; xi++) {
		// Create a Point vector to store the data to give to the data provider
		// for the visualization
		auto myPoints = std::make_shared<std::vector<xolotlViz::Point> >();

		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray + networkSize * xi;
		// Update the concentrations in the network to have physics results
		// (non negative)
		PetscSolver::getNetwork()->updateConcentrationsFromArray(
				gridPointSolution);
		// Get the concentrations from the network
		double concentrations[networkSize];
		double * concentration = &concentrations[0];
		PetscSolver::getNetwork()->fillConcentrationsArray(concentration);

		// A pointer for the clusters used below
		PSICluster * cluster;

		// Loop on Y = V number
		for (int i = 0; i < maxVClusterSize; i++) {
			// Loop on X = He number
			for (int j = 0; j < maxHeVClusterSize - maxVClusterSize; j++) {
				double conc = 0.0;
				// V clusters
				if (j == 0) {
					cluster = (PSICluster *) PetscSolver::getNetwork()->get("V",
							i);
					if (cluster) {
						// Get the ID of the cluster
						int id = cluster->getId() - 1;
						conc = concentration[id];
					}
				}
				// He clusters
				else if (i == 0) {
					cluster = (PSICluster *) PetscSolver::getNetwork()->get(
							"He", j);
					if (cluster) {
						// Get the ID of the cluster
						int id = cluster->getId() - 1;
						conc = concentration[id];
					}
				}
				// HeV clusters
				else {
					cluster =
							(PSICluster *) PetscSolver::getNetwork()->getCompound(
									"HeV", { j, i, 0 });
					if (cluster) {
						// Get the ID of the cluster
						int id = cluster->getId() - 1;
						conc = concentration[id];
					}
				}

				// Create a Point with the concentration as the value
				// and add it to myPoints
				xolotlViz::Point aPoint;
				aPoint.value = conc;
				aPoint.t = time;
				aPoint.x = j;
				aPoint.y = i;
				myPoints->push_back(aPoint);
			}
		}

		// Get the data provider and give it the points
		surfacePlot->getDataProvider()->setPoints(myPoints);
		surfacePlot->getDataProvider()->setDataName("brian");

		// Change the title of the plot
		std::stringstream title;
		title << "Concentration at Depth: " << xi * hx << " nm";
		surfacePlot->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		surfacePlot->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		checkPetscError(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		surfacePlot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "Brian_TS" << timestep << "_D" << xi << ".pnm";
		surfacePlot->write(fileName.str());
	}

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will save 1D plots of one performance timer
 */
PetscErrorCode monitorPerf(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	PetscInt ierr;

	PetscFunctionBeginUser;

	// Get the number of processes
	int cwSize;
    int cwRank;
	MPI_Comm_size(PETSC_COMM_WORLD, &cwSize);
    MPI_Comm_rank(PETSC_COMM_WORLD, &cwRank);

	// Print a warning if only one process
	if (cwSize == 1) {
		std::cout
				<< "You are trying to plot things that don't have any sense!! "
				<< "\nRemove -plot_perf or run in parallel." << std::endl;
		PetscFunctionReturn(0);
	}

    // Obtain the current value of the solve timer.
    //
    // Note that the solve timer keeps a cumultive time,
    // not a per-timestep time.   If you need a per-timestep
    // time, you will want to keep a static or global variable 
    // with the last known timer value, and subtract it from
    // the current timer value each time this monitor function is called.
    //
    // Note also that we restart the timer immediately after sampling
    // its value.  If you feel it is "unfair" to charge the time
    // required for the rank 0 process to produce the output plot
    // against the solve timer, then you should move the start() 
    // call after the code that produces the plot (and probably also
    // put in an MPI_Barrier before starting the timer so that
    // all processes avoid including the time required for rank 0 
    // to produce the plot).  We probably don't want to reset the
    // timer here since the main function is using it to get an
    // overall elapsed time measurement of the solve.
    //
    auto solverTimer = xolotlPerf::getHandlerRegistry()->getTimer("solve");
    solverTimer->stop();
    double solverTimerValue = solverTimer->getValue();
    solverTimer->start();

    // Collect all sampled timer values to rank 0.
    double* allTimerValues = (cwRank == 0) ? new double[cwSize] : NULL;
    MPI_Gather( &solverTimerValue,  // send buffer
                1,                  // number of values to send
                MPI_DOUBLE,         // type of items in send buffer
                allTimerValues,     // receive buffer (only valid at root)
                1,                  // number of values to receive from each process
                MPI_DOUBLE,         // type of items in receive buffer
                0,                  // root of MPI collective operation
                PETSC_COMM_WORLD ); // communicator defining processes involved in the operation

    if( cwRank == 0 )
    {
        auto allPoints = std::make_shared<std::vector<xolotlViz::Point> >();

        for( unsigned int i = 0; i < cwSize; ++i )
        {
            xolotlViz::Point aPoint;
            aPoint.value = allTimerValues[i];
            aPoint.x = cwRank;
            aPoint.t = time;
            allPoints->push_back(aPoint);
        }

        // Provide the data provider the points.
        perfPlot->getDataProvider()->setPoints(allPoints);
		perfPlot->getDataProvider()->setDataName("SolverTimer");

		// Change the title of the plot
		std::ostringstream title;
		title << "Solver timer (s)";
		perfPlot->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::ostringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		perfPlot->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		checkPetscError(ierr);
		// Give the timestep to the label provider
		std::ostringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		perfPlot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::ostringstream fileName;
		fileName << "timer_TS" << timestep << ".pnm";
		perfPlot->write(fileName.str());
    }

    // clean up
    delete[] allTimerValues;

    PetscFunctionReturn(0);
}


PetscErrorCode monitorMaxClusterConc(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int networkSize = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	double *solutionArray, *gridPointSolution, x;
	Vec localSolution;
	int xs, xm, xi;

	PetscFunctionBeginUser;

	// Don't do anything if it was already printed
	if (!printMaxClusterConc) PetscFunctionReturn(0);

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

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Setup some step size variables
	double hx = PetscSolver::getStepSize();

	// Get the maximum size of HeV clusters
	auto psiNetwork = std::dynamic_pointer_cast<PSIClusterReactionNetwork>(
			PetscSolver::getNetwork());
	std::map<std::string, std::string> props = psiNetwork->getProperties();
	int maxHeVClusterSize = std::stoi(props["maxHeVClusterSize"]);
	// Get the maximum size of V clusters
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	// Get the number of He in the max HeV cluster
	int maxHeSize = (maxHeVClusterSize - maxVClusterSize);
	// Get the maximum stable HeV cluster
	PSICluster * maxHeV;
	maxHeV = (PSICluster *) PetscSolver::getNetwork()->getCompound(
			"HeV", { maxHeSize, maxVClusterSize, 0 });

	// Boolean to know if the concentration is too big
	bool maxHeVTooBig = false;

	// Check the concentration of the biggest cluster at each grid point
	for (xi = xs; xi < xs + xm; xi++) {
		// Position
		x = xi * hx;

		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray + networkSize * xi;

		// Update the concentrations in the network
		PetscSolver::getNetwork()->updateConcentrationsFromArray(
				gridPointSolution);

		// Get the concentrations from the network
		double concentrations[networkSize];
		double * concentration = &concentrations[0];
		PetscSolver::getNetwork()->fillConcentrationsArray(
				concentration);
		// Get the concentration of the maximum HeV cluster
		auto maxHeVConc = concentration[maxHeV->getId() - 1];

		if (maxHeVConc > 1.0e-16) maxHeVTooBig = true;
	}

	// Slaves send the information about maxHeVTooBig to the main
	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);
	// Get the current process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Main process
	if (procId == 0) {
		// Loop on the other processes
		for (int i = 1; i < worldSize; i++) {
			int otherBool;
			MPI_Recv(&otherBool, 1, MPI_INT, i, 5, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			// Update maxHeVTooBig
			maxHeVTooBig = maxHeVTooBig || otherBool;
		}

		// Print if maxHeVTooBig is true
		if (maxHeVTooBig) {
			std::cout << std::endl;
			std::cout << "At time step: " << timestep << " and time: " << time
					<< " the biggest cluster: " << *maxHeV
					<< " reached a concentration above 1.0e-16 at at least one grid point."
					<< std::endl << std::endl;

			// Don't print anymore
			printMaxClusterConc = false;
		}

		// Send this information to the other processes
		for (int i = 1; i < worldSize; i++) {
			int printMaxClusterConcInt = (int) printMaxClusterConc;
			MPI_Send(&printMaxClusterConcInt, 1, MPI_INT, i, 5, MPI_COMM_WORLD);
		}
	}
	// Other processes
	else {
		// Send the maxHeVTooBig value
		int maxHeVTooBigInt = (int) maxHeVTooBig;
		MPI_Send(&maxHeVTooBigInt, 1, MPI_INT, 0, 5, MPI_COMM_WORLD);

		// Receive the printMaxClusterConc value
		MPI_Recv(&printMaxClusterConc, 1, MPI_INT, 0, 5, MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);
	}

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will compute the flux of interstitials
 * at the surface
 */
PetscErrorCode monitorInterstitial(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int networkSize = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	double *solutionArray, *gridPointSolution, x;
	Vec localSolution;
	int xs, xm, xi;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

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

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Return if the grid point 1 is not on this process
	if (xs > 1) PetscFunctionReturn(0);

	// Setup some step size variables
	double hx = PetscSolver::getStepSize();
	double sx = 1.0 / (hx * hx);

	// Get the concentrations at xi = 1
	xi = 1;
	gridPointSolution = solutionArray + networkSize * xi;
	// Update the concentrations in the network to have physics results
	// (non negative)
	PetscSolver::getNetwork()->updateConcentrationsFromArray(
			gridPointSolution);
	// Get the concentrations from the network
	double concentrations[networkSize];
	double * concentration = &concentrations[0];
	PetscSolver::getNetwork()->fillConcentrationsArray(concentration);

	// Get the delta time from the previous timestep to this timestep
	double dt = time - previousTime;

	// Compute the total density of intersitials that escaped from the
	// surface since last timestep using the stored flux
	double nInterstitial = previousFlux * dt;

//	// Uncomment to write the interstitial flux in a file
//	std::ofstream outputFile;
//	outputFile.open("interstitialOut.txt", ios::app);
//	outputFile << time << " "
//			<< nInterstitial << std::endl;
//	outputFile.close();

	// Initialize the value for the flux
	double newFlux = 0.0;

	// Get all the interstitial clusters
	auto interstitials = PetscSolver::getNetwork()->getAll("I");
	// Loop on them
	for (int i = 0; i < interstitials.size(); i++) {
		// Get the cluster
		auto cluster = (PSICluster *) interstitials.at(i);
		// Get its id and concentration
		int id = cluster->getId() - 1;
		double conc = concentration[id];
		// Get its size and diffusion coefficient
		int size = cluster->getSize();
		double coef = cluster->getDiffusionCoefficient();

		// Compute the flux
		newFlux += (double) size * sx * coef * conc;
	}

	previousFlux = newFlux;

    PetscFunctionReturn(0);
}

/**
 * This operation sets up a monitor that will call monitorSolve
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
PetscErrorCode setupPetscMonitor(TS ts) {
	PetscErrorCode ierr;

	//! The xolotlViz handler registry
	auto vizHandlerRegistry = xolotlViz::getVizHandlerRegistry();

	// Flags to launch the monitors or not
	PetscBool flag2DPlot, flag1DPlot, flagSeries, flagPerf, flagRetention,
			flagStatus, flagMaxClusterConc, flagInterstitial;

	// Check the option -plot_perf
	ierr = PetscOptionsHasName(NULL, "-plot_perf", &flagPerf);
	checkPetscError(ierr);

	// Check the option -plot_series
	ierr = PetscOptionsHasName(NULL, "-plot_series", &flagSeries);
	checkPetscError(ierr);

	// Check the option -plot_1d
	ierr = PetscOptionsHasName(NULL, "-plot_1d", &flag1DPlot);
	checkPetscError(ierr);

	// Check the option -plot_2d
	ierr = PetscOptionsHasName(NULL, "-plot_2d", &flag2DPlot);
	checkPetscError(ierr);

	// Check the option -helium_retention
	ierr = PetscOptionsHasName(NULL, "-helium_retention", &flagRetention);
	checkPetscError(ierr);

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, "-start_stop", &flagStatus);
	checkPetscError(ierr);

	// Check the option -max_cluster_conc
	ierr = PetscOptionsHasName(NULL, "-max_cluster_conc", &flagMaxClusterConc);
	checkPetscError(ierr);

	// Check the option -interstitial_diff
	ierr = PetscOptionsHasName(NULL, "-interstitial_diff", &flagInterstitial);
	checkPetscError(ierr);

	// Set the monitor to save 1D plot of one concentration
	if (flag1DPlot) {
		// Create a ScatterPlot
		plot = vizHandlerRegistry->getPlot("scatterPlot",
				xolotlViz::PlotType::SCATTER);

		// Create and set the label provider
		auto labelProvider = std::make_shared<xolotlViz::LabelProvider>(
				"labelProvider");
		labelProvider->axis1Label = "x Position on the Grid";
		labelProvider->axis2Label = "Concentration";

		// Give it to the plot
		plot->setLabelProvider(labelProvider);

		// Create the data provider
		auto dataProvider = std::make_shared<xolotlViz::CvsXDataProvider>(
				"dataProvider");

		// Give it to the plot
		plot->setDataProvider(dataProvider);

		// monitorSolve will be called at each timestep
		ierr = TSMonitorSet(ts, monitorScatter, NULL, NULL);
		checkPetscError(ierr);
	}

	// Set the monitor to save 1D plot of many concentrations
	if (flagSeries) {
		// Create a ScatterPlot
		seriesPlot = vizHandlerRegistry->getPlot("seriesPlot",
				xolotlViz::PlotType::SERIES);

		// set the log scale
		seriesPlot->setLogScale();

		// Create and set the label provider
		auto labelProvider = std::make_shared<xolotlViz::LabelProvider>(
				"labelProvider");
		labelProvider->axis1Label = "x Position on the Grid";
		labelProvider->axis2Label = "Concentration";

		// Give it to the plot
		seriesPlot->setLabelProvider(labelProvider);

		// Network size
		const int networkSize = PetscSolver::getNetwork()->size();

		// To plot a maximum of 18 clusters of the whole benchmark
		const int loopSize = std::min(18, networkSize);

		// Create a data provider for each cluster in the network
		for (int i = 0; i < loopSize; i++) {
			// Set the name for Identifiable
			std::stringstream dataProviderName;
			dataProviderName << "dataprovider" << i;
			// Create the data provider
			auto dataProvider = std::make_shared<xolotlViz::CvsXDataProvider>(
					dataProviderName.str());

			// Give it to the plot
			seriesPlot->addDataProvider(dataProvider);
		}

		// monitorSolve will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSeries, NULL, NULL);
		checkPetscError(ierr);
	}

	// Set the monitor to save surface plots of clusters concentration
	// for each depth
	if (flag2DPlot) {
		// Create a SurfacePlot
		surfacePlot = vizHandlerRegistry->getPlot("surfacePlot",
				xolotlViz::PlotType::SURFACE);

		// Create and set the label provider
		auto labelProvider = std::make_shared<xolotlViz::LabelProvider>(
				"labelProvider");
		labelProvider->axis1Label = "He number";
		labelProvider->axis2Label = "V number";
		labelProvider->axis3Label = "Concentration";

		// Give it to the plot
		surfacePlot->setLabelProvider(labelProvider);

		// Create the data provider
		auto dataProvider = std::make_shared<xolotlViz::CvsXYDataProvider>(
				"dataProvider");

		// Give it to the plot
		surfacePlot->setDataProvider(dataProvider);

		// monitorSeries will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSurface, NULL, NULL);
		checkPetscError(ierr);
	}

	// Set the monitor to save performance plots (has to be in parallel)
	if (flagPerf) {
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

		// monitorPerf will be called at each timestep
		ierr = TSMonitorSet(ts, monitorPerf, NULL, NULL);
		checkPetscError(ierr);

	}

	// Set the monitor to compute the helium fluence for the retention calculation
	if (flagRetention) {
		// Get all the helium clusters
		auto heClusters = PetscSolver::getNetwork()->getAll(heType);

		// Get all the helium-vacancy clusters
		auto heVClusters = PetscSolver::getNetwork()->getAll(heVType);

		// Loop on the helium clusters
		for (int i = 0; i < heClusters.size(); i++) {
			auto cluster = (PSICluster *) heClusters[i];
			int id = cluster->getId() - 1;
			// Add the Id to the vector
			indices.push_back(id);
			// Add the number of heliums of this cluster to the weight
			weight.push_back(cluster->getSize());
		}

		// Loop on the helium-vacancy clusters
		for (int i = 0; i < heVClusters.size(); i++) {
			auto cluster = (PSICluster *) heVClusters[i];
			int id = cluster->getId() - 1;
			// Add the Id to the vector
			indices.push_back(id);
			// Add the number of heliums of this cluster to the weight
			auto comp = cluster->getComposition();
			weight.push_back(comp[heType]);
		}

		if (indices.size() == 0) {
			throw std::string(
					"PetscSolver Exception: Cannot compute the retention because there is no helium or helium-vacancy cluster in the network.");
		}

		// computeHeliumFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumFluence, NULL, NULL);
		checkPetscError(ierr);

		// computeHeliumRetention will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumRetention, NULL, NULL);
		checkPetscError(ierr);

//		// Uncomment to clear the file where the retention will be written
//		std::ofstream outputFile;
//		outputFile.open("retentionOut.txt");
//		outputFile.close();
	}

	// Set only the monitor to compute the fluence in the case the retention
	// option is not used but the fluence one is
	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = PetscSolver::getFluxHandler();
	// Get the fluence option
	bool heFluenceOption = fluxHandler->getUsingMaxHeFluence();
	if (!flagRetention && heFluenceOption) {
		// computeHeliumFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumFluence, NULL, NULL);
		checkPetscError(ierr);
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr = PetscOptionsGetInt(NULL, "-start_stop", &stride, &flag);
		checkPetscError(ierr);
		if (!flag)
			stride = 1;

		// Network size
		const int networkSize = PetscSolver::getNetwork()->size();
		PetscInt Mx;
		PetscErrorCode ierr;

		// Get the da from ts
		DM da;
		ierr = TSGetDM(ts, &da);
		checkPetscError(ierr);

		// Get the size of the total grid
		ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE);
		checkPetscError(ierr);

		// Initialize the HDF5 file for all the processes
		xolotlCore::HDF5Utils::initializeFile(outputFileName, networkSize, Mx);

		// Setup step size variable
		double hx = PetscSolver::getStepSize();

		// Get the physical dimension of the grid
		int dimension = (Mx - 1) * hx;

		// Get the refinement of the grid
		PetscInt refinement = 0;
		ierr = PetscOptionsGetInt(NULL, "-da_refine", &refinement, &flag);
		checkPetscError(ierr);
		if (!flag)
			refinement = 0;

		// Save the header in the HDF5 file
		xolotlCore::HDF5Utils::fillHeader(dimension, refinement);

		// Save the network in the HDF5 file
		xolotlCore::HDF5Utils::fillNetwork(PetscSolver::getNetwork());

		// Finalize the HDF5 file
		xolotlCore::HDF5Utils::finalizeFile();

		// startStop will be called at each timestep
		ierr = TSMonitorSet(ts, startStop, NULL, NULL);
		checkPetscError(ierr);
	}

	// Set the monitor to output information about when the maximum stable HeV
	// cluster in the network first becomes greater than 1.0e-16
	if (flagMaxClusterConc) {
		// monitorMaxClusterConc will be called at each timestep
		ierr = TSMonitorSet(ts, monitorMaxClusterConc, NULL, NULL);
		checkPetscError(ierr);
	}

	// Set the monitor on the outgoing flux of interstitials at the surface
	if (flagInterstitial) {
		// monitorInterstitial will be called at each timestep
		ierr = TSMonitorSet(ts, monitorInterstitial, NULL, NULL);
		checkPetscError(ierr);

//		// Uncomment to clear the file where the interstitial flux will be written
//		std::ofstream outputFile;
//		outputFile.open("interstitialOut.txt");
//		outputFile.close();
	}

	// Set the monitor to simply change the previous time to the new time
	if (flagRetention || heFluenceOption || flagInterstitial) {
		// monitorTime will be called at each timestep
		ierr = TSMonitorSet(ts, monitorTime, NULL, NULL);
		checkPetscError(ierr);
	}

	PetscFunctionReturn(0);
}

}

/* end namespace xolotlSolver */
