// Includes
#include "PetscSolver.h"
#include <HandlerRegistryFactory.h>
#include <VizHandlerRegistryFactory.h>
#include <PlotType.h>
#include <CvsXDataProvider.h>
#include <CvsXYDataProvider.h>
#include <LabelProvider.h>
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
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

//! Physical length of the grid
int xGridLength = 8;

//! The double that will store the accumulation of helium flux.
double heliumFluence = 0.0;

//! How often HDF5 file is written
PetscInt stride = 0;

//! HDF5 output file name
std::string outputFileName = "xolotlStop.h5";

/**
 * This is a monitoring method that will save an hdf5 file at each time step.
 * HDF5 is handling the parallel part, so no call to MPI here.
 */
static PetscErrorCode startStop(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int networkSize = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	PetscReal *solutionArray, *gridPointSolution;
	Vec localSolution;
	PetscInt xs, xm, Mx;

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
	double hx = (double) xGridLength / (PetscReal)(Mx - 1);

	// Open the already created HDF5 file
	xolotlCore::HDF5Utils::openFile(outputFileName);

	// Get the physical dimension of the grid
	int dimension = (Mx - 1) * hx;

	// Get the current time step
	PetscReal currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	checkPetscError(ierr);

	// Add a concentration sub group
	xolotlCore::HDF5Utils::addConcentrationSubGroup(timestep, networkSize,
			time, currentTimeStep);

	// Loop on the full grid
	for (int xi = 0; xi < Mx; xi++) {
		// Size of the concentration that will be stored
		int concSize = -1;
		// Vector for the concentrations
		std::vector< std::vector<double> > concVector;

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
				if (i == procId) continue;

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
		if (concSize == 0) continue;

		// All processes must create the dataset
		xolotlCore::HDF5Utils::addConcentrationDataset(xi, concSize);

		MPI_Barrier(MPI_COMM_WORLD);

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
static PetscErrorCode heliumRetention(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int size = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	PetscInt xs, xm, Mx;

	PetscFunctionBeginUser;

	// Variable to represent the real, or current, time
	PetscReal realTime;
	ierr = TSGetTime(ts, &realTime);

	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = PetscSolver::getFluxHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
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
	double hx = (double) xGridLength / (PetscReal)(Mx - 1);

	// Get the helium cluster
	auto heCluster = (PSICluster *) PetscSolver::getNetwork()->get("He", 1);

	// Exit if there is no helium cluster in the network
	if (!heCluster)
		PetscFunctionReturn(0);

	// Get the composition of the cluster
	auto thisComp = heCluster->getComposition();
	std::vector<int> compVec = { thisComp["He"], thisComp["V"], thisComp["I"] };

	// Define the reactant ID
	int reactantIndex = heCluster->getId() - 1;

	// Loop on the grid
	for (int xi = xs; xi < xs + xm; xi++) {
		// Actual position in nm
		double x = xi * hx;

		// Vector representing the position at which the flux will be calculated
		// Currently we are only in 1D
		std::vector<double> gridPosition = { x, 0, 0 };

		// Calculate the incident flux
		auto incidentFlux = fluxHandler->getIncidentFlux(compVec, gridPosition,
				realTime);

		// And add it to the fluence
		heliumFluence += 10000.0 * std::max(0.0, incidentFlux) * time;
	}

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will save 1D plots of one concentration
 */
static PetscErrorCode monitorScatter(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int networkSize = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	PetscReal *solutionArray, *gridPointSolution, x, hx;
	Vec localSolution;
	PetscInt xs, xm, Mx;
	int xi, i;

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
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr);
	// Setup some step size variables
	hx = (double) xGridLength / (PetscReal)(Mx - 1);

	// Choice of the cluster to be plotted
	int iCluster = 7;

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
			MPI_Recv(&localSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			// Loop on their grid
			for (int k = 0; k < localSize; k++) {
				// Get the position
				MPI_Recv(&x, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

				// and the concentration
				double conc = 0.0;
				MPI_Recv(&conc, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
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
		MPI_Send(&xm, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

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
			MPI_Send(&x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

			// Send the value of the concentration to the master process
			MPI_Send(&concentration[iCluster], 1, MPI_DOUBLE, 0, 0,
					MPI_COMM_WORLD);
		}
	}

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will save 1D plots of many concentrations
 */
static PetscErrorCode monitorSeries(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int networkSize = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	PetscReal *solutionArray, *gridPointSolution, x, hx;
	Vec localSolution;
	PetscInt xs, xm, Mx;
	int xi, i;

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
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr);
	// Setup some step size variables
	hx = (double) xGridLength / (PetscReal)(Mx - 1);

	// To plot a maximum of 18 clusters of the whole benchmark
	const int loopSize = std::min(18, networkSize);

	if (procId == 0) {
		// Create a Point vector to store the data to give to the data provider
		// for the visualization
		std::vector < std::vector<xolotlViz::Point> > myPoints(loopSize);

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
			MPI_Recv(&localSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			// Loop on their grid
			for (int k = 0; k < localSize; k++) {
				// Get the position
				MPI_Recv(&x, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

				for (int j = 0; j < loopSize; j++) {
					// and the concentrations
					double conc;
					MPI_Recv(&conc, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
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
			auto thePoints = std::make_shared < std::vector<xolotlViz::Point>
					> (myPoints[i]);
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
		MPI_Send(&xm, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

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
			MPI_Send(&x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

			for (int i = 0; i < loopSize; i++) {
				// Send the value of the concentrations to the master process
				MPI_Send(&concentration[i], 1, MPI_DOUBLE, 0, 0,
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
static PetscErrorCode monitorSurface(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int networkSize = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	PetscReal *solutionArray, *gridPointSolution, x, hx;
	Vec localSolution;
	PetscInt xs, xm, Mx;
	int xi, i;

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
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr);
	// Setup some step size variables
	hx = (double) xGridLength / (PetscReal)(Mx - 1);

	// Get the maximum size of HeV clusters
	auto psiNetwork = std::dynamic_pointer_cast < PSIClusterReactionNetwork
			> (PetscSolver::getNetwork());
	std::map < std::string, std::string > props = psiNetwork->getProperties();
	int maxHeVClusterSize = std::stoi(props["maxHeVClusterSize"]);

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
		for (int i = 0; i < maxHeVClusterSize; i++) {
			// Loop on X = He number
			for (int j = 0; j < maxHeVClusterSize; j++) {
				double conc = 0.0;
				// V clusters
				if (j == 0) {
					cluster = (PSICluster *) PetscSolver::getNetwork()->get("V", i);
					if (cluster) {
						// Get the ID of the cluster
						int id = cluster->getId() - 1;
						conc = concentration[id];
					}
				}
				// He clusters
				else if (i == 0) {
					cluster = (PSICluster *) PetscSolver::getNetwork()->get("He", j);
					if (cluster) {
						// Get the ID of the cluster
						int id = cluster->getId() - 1;
						conc = concentration[id];
					}
				}
				// HeV clusters
				else {
					cluster = (PSICluster *) PetscSolver::getNetwork()->getCompound("HeV", {
									j, i, 0 });
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
static PetscErrorCode monitorPerf(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	PetscInt ierr;

	PetscFunctionBeginUser;

	// Get the number of processes
	int size;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);

	// Print a warning if only one process
	if (size == 1) {
		std::cout
				<< "You are trying to plot things that don't have any sense!! "
				<< "\nRemove -plot_perf or run in parallel." << std::endl;
		PetscFunctionReturn(0);
	}

	// Get the current process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the solve timer
	auto solverTimer = xolotlPerf::getHandlerRegistry()->getTimer("solve");

	// Stop it to access its value
	solverTimer->stop();

	// Master process
	if (procId == 0) {

		// Create a Point vector to store the data to give to the data provider
		// for the visualization
		auto myPoints = std::make_shared<std::vector<xolotlViz::Point> >();

		// Give it the value for procId = 0
		xolotlViz::Point aPoint;
		aPoint.value = solverTimer->getValue();
		aPoint.t = time;
		aPoint.x = procId;
		myPoints->push_back(aPoint);

		// Loop on all the other processes
		for (int i = 1; i < size; i++) {
			double counter = 0.0;

			// Receive the value from the other processes
			MPI_Recv(&counter, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			// Give it the value for procId = i
			aPoint.value = counter;
			aPoint.t = time;
			aPoint.x = i;
			myPoints->push_back(aPoint);
		}

		// Get the data provider and give it the points
		perfPlot->getDataProvider()->setPoints(myPoints);

		perfPlot->getDataProvider()->setDataName("SolverTimer");

		// Change the title of the plot
		std::stringstream title;
		title << "Solver timer (s)";
		perfPlot->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		perfPlot->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		checkPetscError(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		perfPlot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "timer_TS" << timestep << ".pnm";
		perfPlot->write(fileName.str());
	}

	else {
		double counter = solverTimer->getValue();

		// Send the value of the timer to the master process
		MPI_Send(&counter, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	// Restart the timer
	solverTimer->start();

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

	// Get the physical length of the grid
	PetscBool flg;
	PetscInt length;
	PetscOptionsGetInt(NULL, "-da_grid_x", &length, &flg);
	if (flg) xGridLength = (int) length;
	else xGridLength = 8;

	// Flags to launch the monitors or not
	PetscBool flag2DPlot, flag1DPlot, flagSeries, flagPerf, flagRetention,
			flagStatus;

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

	// Use this instead of trying every flag
	int bigFlag = flag1DPlot + flag2DPlot + flagSeries + flagPerf
			+ flagRetention + flagStatus;

	// Don't do anything if no option is set
	if (!bigFlag)
		PetscFunctionReturn(0);

	// Set the monitor to save 1D plot of one concentration
	if (flag1DPlot) {
		// Create a ScatterPlot
		plot = vizHandlerRegistry->getPlot("scatterPlot",
				xolotlViz::PlotType::SCATTER);

		// Create and set the label provider
		auto labelProvider = std::make_shared < xolotlViz::LabelProvider
				> ("labelProvider");
		labelProvider->axis1Label = "x Position on the Grid";
		labelProvider->axis2Label = "Concentration";

		// Give it to the plot
		plot->setLabelProvider(labelProvider);

		// Create the data provider
		auto dataProvider = std::make_shared < xolotlViz::CvsXDataProvider
				> ("dataProvider");

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
		auto labelProvider = std::make_shared < xolotlViz::LabelProvider
				> ("labelProvider");
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
			auto dataProvider = std::make_shared < xolotlViz::CvsXDataProvider
					> (dataProviderName.str());

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
		auto labelProvider = std::make_shared < xolotlViz::LabelProvider
				> ("labelProvider");
		labelProvider->axis1Label = "He number";
		labelProvider->axis2Label = "V number";
		labelProvider->axis3Label = "Concentration";

		// Give it to the plot
		surfacePlot->setLabelProvider(labelProvider);

		// Create the data provider
		auto dataProvider = std::make_shared < xolotlViz::CvsXYDataProvider
				> ("dataProvider");

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
		auto labelProvider = std::make_shared < xolotlViz::LabelProvider
				> ("labelProvider");
		labelProvider->axis1Label = "Process ID";
		labelProvider->axis2Label = "Solver Time";

		// Give it to the plot
		perfPlot->setLabelProvider(labelProvider);

		// Create the data provider
		auto dataProvider = std::make_shared < xolotlViz::CvsXDataProvider
				> ("dataProvider");

		// Give it to the plot
		perfPlot->setDataProvider(dataProvider);

		// monitorPerf will be called at each timestep
		ierr = TSMonitorSet(ts, monitorPerf, NULL, NULL);
		checkPetscError(ierr);

	}

	// Set the monitor to compute the helium fluence for the retention calculation
	if (flagRetention) {
		// heliumRetention will be called at each timestep
		ierr = TSMonitorSet(ts, heliumRetention, NULL, NULL);
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
		double hx = (double) xGridLength / (PetscReal) (Mx - 1);

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

	PetscFunctionReturn(0);
}

/**
 * This operation computes and prints the helium retention.
 * @param ts The time stepper.
 * @param C The vector of solution.
 */
void computeRetention(TS ts, Vec C) {
	PetscErrorCode ierr;
	// Network size
	const int size = PetscSolver::getNetwork()->size();

	// Declare the vector that will store the Id of the helium clusters
	std::vector<int> indices;

	// Declare the vector that will store the weight of the helium clusters
	// (their He composition)
	std::vector<int> weight;

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
		return;
	}

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the array of concentration
	PetscReal *solutionArray;
	ierr = DMDAVecGetArray(da, C, &solutionArray);
	checkPetscError(ierr);

	//Get local grid boundaries
	PetscInt xs, xm;
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Get the size of the total grid
	PetscInt Mx;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Setup step size variable
	double hx = (double) xGridLength / (PetscReal) (Mx - 1);

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
			double otherFluence = 0.0;

			// Receive the value from the other processes
			MPI_Recv(&otherConcentration, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);
			MPI_Recv(&otherFluence, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			// Add them to the master one
			heConcentration += otherConcentration;
			heliumFluence += otherFluence;
		}

		// Get the final time
		PetscReal time;
		ierr = TSGetTime(ts, &time);
		checkPetscError(ierr);

		// Print the result
		std::cout << "Final time: " << time << std::endl;
		std::cout << "Helium retention = "
				<< 100.0 * (heConcentration / heliumFluence) << " %"
				<< std::endl;
		std::cout << "Helium concentration = " << heConcentration << std::endl;
		std::cout << "Helium fluence = " << heliumFluence << std::endl;
	}

	else {
		// Send the value of the timer to the master process
		MPI_Send(&heConcentration, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&heliumFluence, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	return;
}

}

/* end namespace xolotlSolver */
