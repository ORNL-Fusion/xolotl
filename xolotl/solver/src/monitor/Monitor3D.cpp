// Includes
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <xolotl/solver/PetscSolver.h>
#include <xolotl/perf/xolotlPerf.h>
#include <xolotl/factory/viz/VizHandlerRegistryFactory.h>
#include <xolotl/viz/PlotType.h>
#include <xolotl/viz/dataprovider/CvsXDataProvider.h>
#include <xolotl/viz/dataprovider/CvsXYDataProvider.h>
#include <xolotl/viz/LabelProvider.h>
#include <xolotl/core/Constants.h>
#include <xolotl/util/MathUtils.h>
#include <xolotl/util/RandomNumberGenerator.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/io/XFile.h>
#include <xolotl/solver/monitor/Monitor.h>

namespace xolotl {
namespace solver {
namespace monitor {

// Declaration of the functions defined in Monitor.cpp
extern PetscErrorCode checkTimeStep(TS ts);
extern PetscErrorCode monitorTime(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode computeFluence(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode monitorPerf(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);

// Declaration of the variables defined in Monitor.cpp
extern std::shared_ptr<viz::IPlot> perfPlot;
extern double previousTime;
extern double timeStepThreshold;

//! How often HDF5 file is written
PetscReal hdf5Stride3D = 0.0;
//! Previous time for HDF5
PetscInt hdf5Previous3D = 0;
//! HDF5 output file name
std::string hdf5OutputName3D = "xolotlStop.h5";
//! The pointer to the 2D plot used in MonitorSurfaceXY3D.
std::shared_ptr<viz::IPlot> surfacePlotXY3D;
//! The pointer to the 2D plot used in MonitorSurfaceXZ3D.
std::shared_ptr<viz::IPlot> surfacePlotXZ3D;
//! The variable to store the interstitial flux at the previous time step.
std::vector<std::vector<double> > previousIFlux3D;
//! The variable to store the total number of interstitials going through the surface.
std::vector<std::vector<double> > nInterstitial3D;
//! The variable to store the sputtering yield at the surface.
double sputteringYield3D = 0.0;
// The vector of depths at which bursting happens
std::vector<std::tuple<int, int, int> > depthPositions3D;
// The vector of ids for diffusing interstitial clusters
std::vector<int> iClusterIds3D;

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "startStop3D")
/**
 * This is a monitoring method that will update an hdf5 file at each time step.
 */
PetscErrorCode startStop3D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void*) {
	// Initial declarations
	PetscErrorCode ierr;
	const double ****solutionArray, *gridPointSolution;
	PetscInt xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Compute the dt
	double dt = time - previousTime;

	// Don't do anything if it is not on the stride
	if (((int) ((time + dt / 10.0) / hdf5Stride3D) <= hdf5Previous3D)
			&& timestep > 0)
		PetscFunctionReturn(0);

	// Update the previous time
	if ((int) ((time + dt / 10.0) / hdf5Stride3D) > hdf5Previous3D)
		hdf5Previous3D++;

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
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network and dof
	auto &network = solverHandler.getNetwork();
	const int dof = network.getDOF();

	// Create an array for the concentration
	double concArray[dof][2];

	// Get the vector of positions of the surface
	std::vector<std::vector<int> > surfaceIndices;
	for (PetscInt i = 0; i < My; i++) {
		// Create a temporary vector
		std::vector<int> temp;
		for (PetscInt j = 0; j < Mz; j++) {
			temp.push_back(solverHandler.getSurfacePosition(i, j));
		}
		// Add the temporary vector to the vector of surface indices
		surfaceIndices.push_back(temp);
	}

	// Open the existing HDF5 file.
	io::XFile checkpointFile(hdf5OutputName3D, PETSC_COMM_WORLD,
			io::XFile::AccessMode::OpenReadWrite);

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	CHKERRQ(ierr);

	// Add a concentration sub group
	auto concGroup = checkpointFile.getGroup<
			io::XFile::ConcentrationGroup>();
	assert(concGroup);
	auto tsGroup = concGroup->addTimestepGroup(timestep, time, previousTime,
			currentTimeStep);

	if (solverHandler.moveSurface()) {
		// Write the surface positions in the concentration sub group
		tsGroup->writeSurface3D(surfaceIndices, nInterstitial3D,
				previousIFlux3D);
	}

	// Loop on the full grid
	for (PetscInt k = 0; k < Mz; k++) {
		for (PetscInt j = 0; j < My; j++) {
			for (PetscInt i = 0; i < Mx; i++) {
				// Wait for all the processes
				MPI_Barrier(PETSC_COMM_WORLD);

				// Size of the concentration that will be stored
				int concSize = -1;
				// To save which proc has the information
				int concId = 0;
				// To know which process should write
				bool write = false;

				// If it is the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym && k >= zs
						&& k < zs + zm) {
					write = true;
					// Get the pointer to the beginning of the solution data for this grid point
					gridPointSolution = solutionArray[k][j][i];

					// Loop on the concentrations
					for (int l = 0; l < dof + 1; l++) {
						if (std::fabs(gridPointSolution[l]) > 1.0e-16) {
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

				// All processes create the dataset and fill it
				tsGroup->writeConcentrationDataset(concSize, concArray, write,
						i, j, k);
			}
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeHeliumRetention3D")
/**
 * This is a monitoring method that will compute the helium retention.
 */
PetscErrorCode computeHeliumRetention3D(TS ts, PetscInt, PetscReal time,
		Vec solution, void*) {
	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm, ys, ym, zs, zm, Mx, My, Mz;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = solverHandler.getFluxHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid in the x direction
	auto grid = solverHandler.getXGrid();

	// Get the network
	using NetworkType =
	core::network::PSIReactionNetwork<core::network::PSIFullSpeciesList>;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	auto &network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	const int dof = network.getDOF();

	// Setup step size variables
	double hy = solverHandler.getStepSizeY();
	double hz = solverHandler.getStepSizeZ();

	// Get the array of concentration
	double ****solutionArray, *gridPointSolution;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration over the grid
	double heConcentration = 0.0, dConcentration = 0.0, tConcentration = 0.0;

	// Loop on the grid
	for (PetscInt zk = zs; zk < zs + zm; zk++) {
		for (PetscInt yj = ys; yj < ys + ym; yj++) {

			// Get the surface position
			int surfacePos = solverHandler.getSurfacePosition(yj, zk);

			for (PetscInt xi = xs; xi < xs + xm; xi++) {

				// Boundary conditions
				if (xi < surfacePos + solverHandler.getLeftOffset()
						|| xi >= Mx - solverHandler.getRightOffset())
					continue;

				// Get the pointer to the beginning of the solution data for
				// this grid point
				gridPointSolution = solutionArray[zk][yj][xi];

				double hx = grid[xi + 1] - grid[xi];

				using HostUnmanaged =
				Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
				auto hConcs = HostUnmanaged(gridPointSolution, dof);
				auto dConcs = Kokkos::View<double*>("Concentrations", dof);
				deep_copy(dConcs, hConcs);

				// Get the total concentrations at this grid point
				heConcentration += network.getTotalAtomConcentration(dConcs,
						Spec::He, 0) * hx * hy * hz;
				dConcentration += network.getTotalAtomConcentration(dConcs,
						Spec::D, 0) * hx * hy * hz;
				tConcentration += network.getTotalAtomConcentration(dConcs,
						Spec::T, 0) * hx * hy * hz;
			}
		}
	}

	// Get the current process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Determine total concentrations for He, D, T.
	std::array<double, 3> myConcData { heConcentration, dConcentration,
			tConcentration };
	std::array<double, 3> totalConcData;

	MPI_Reduce(myConcData.data(), totalConcData.data(), myConcData.size(),
	MPI_DOUBLE,
	MPI_SUM, 0, PETSC_COMM_WORLD);

	// Extract total He, D, T concentrations.  Values are valid only on rank 0.
	double totalHeConcentration = totalConcData[0];
	double totalDConcentration = totalConcData[1];
	double totalTConcentration = totalConcData[2];

	// Master process
	if (procId == 0) {

		// Get the total size of the grid rescale the concentrations
		PetscInt Mx, My, Mz;
		ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
		CHKERRQ(ierr);

		// Compute the total surface irradiated by the helium flux
		double surface = (double) (My * Mz) * hy * hz;

		// Rescale the concentration
		totalHeConcentration = totalHeConcentration / surface;
		totalDConcentration = totalDConcentration / surface;
		totalTConcentration = totalTConcentration / surface;

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
		outputFile.open("retentionOut.txt", std::ios::app);
		outputFile << fluence << " " << totalHeConcentration << " "
				<< totalDConcentration << " " << totalTConcentration
				<< std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeXenonRetention3D")
/**
 * This is a monitoring method that will compute the xenon retention
 */
PetscErrorCode computeXenonRetention3D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void*) {

	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm, ys, ym, zs, zm;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	CHKERRQ(ierr);

	// Get the total size of the grid
	PetscInt Mx, My, Mz;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// Setup step size variables
	double hy = solverHandler.getStepSizeY();
	double hz = solverHandler.getStepSizeZ();

	using NetworkType =
	core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Degrees of freedom is the total number of clusters in the network
	auto &network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	const int dof = network.getDOF();

	// Get the array of concentration
	PetscReal ****solutionArray, *gridPointSolution;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration and other values over the grid
	double xeConcentration = 0.0, bubbleConcentration = 0.0, radii = 0.0,
			partialBubbleConcentration = 0.0, partialRadii = 0.0;

	// Get the minimum size for the radius
	auto minSizes = solverHandler.getMinSizes();

	// Loop on the grid
	for (PetscInt zk = zs; zk < zs + zm; zk++) {
		for (PetscInt yj = ys; yj < ys + ym; yj++) {
			for (PetscInt xi = xs; xi < xs + xm; xi++) {

				// Get the pointer to the beginning of the solution data for this grid point
				gridPointSolution = solutionArray[zk][yj][xi];

				using HostUnmanaged =
				Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
				auto hConcs = HostUnmanaged(gridPointSolution, dof);
				auto dConcs = Kokkos::View<double*>("Concentrations", dof);
				deep_copy(dConcs, hConcs);

				double hx = grid[xi + 1] - grid[xi];

				// Get the concentrations
				xeConcentration += network.getTotalAtomConcentration(dConcs,
						Spec::Xe, 0) * hx * hy * hz;
				bubbleConcentration += network.getTotalConcentration(dConcs,
						Spec::Xe, 0) * hx * hy * hz;
				radii += network.getTotalRadiusConcentration(dConcs, Spec::Xe,
						0) * hx * hy * hz;
				partialBubbleConcentration += network.getTotalConcentration(
						dConcs, Spec::Xe, minSizes[0]) * hx * hy * hz;
				partialRadii += network.getTotalRadiusConcentration(dConcs,
						Spec::Xe, minSizes[0]) * hx * hy * hz;
			}
		}
	}

	// Get the current process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Sum all the concentrations through MPI reduce
	std::array<double, 5> myConcData { xeConcentration, bubbleConcentration,
			radii, partialBubbleConcentration, partialRadii };
	std::array<double, 5> totalConcData { 0.0, 0.0, 0.0, 0.0, 0.0 };
	MPI_Reduce(myConcData.data(), totalConcData.data(), myConcData.size(),
	MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);

	// Master process
	if (procId == 0) {
		// Compute the total surface irradiated
		double surface = (double) My * hy * (double) Mz * hz;

		totalConcData[0] = totalConcData[0] / surface;

		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		std::cout << "Xenon concentration = " << totalConcData[0] << std::endl
				<< std::endl;

		// Make sure the average partial radius makes sense
		double averagePartialRadius = 0.0;
		if (totalConcData[3] > 1.e-16) {
			averagePartialRadius = totalConcData[4] / totalConcData[3];
		}

		// Uncomment to write the retention and the fluence in a file
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt", std::ios::app);
		outputFile << time << " " << totalConcData[0] << " "
				<< totalConcData[2] / totalConcData[1] << " "
				<< averagePartialRadius << std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorSurfaceXY3D")
/**
 * This is a monitoring method that will save 2D plots of the concentration of
 * a specific cluster at each grid point on the XY surface, integrating over Z.
 */
PetscErrorCode monitorSurfaceXY3D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void*) {
	// Initial declarations
	PetscErrorCode ierr;
	const double ****solutionArray, *gridPointSolution;
	PetscInt xs, xm, Mx, ys, ym, My, zs, zm, Mz;
	double x, y;

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
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the physical grid in the x direction
	auto grid = solverHandler.getXGrid();

	// Setup step size variables
	double hy = solverHandler.getStepSizeY();

	// Choice of the cluster to be plotted
	int iCluster = 0;

	// Create a Point vector to store the data to give to the data provider
	// for the visualization
	auto myPoints = std::make_shared<std::vector<viz::dataprovider::Point> >();
	// Create a point here so that it is not created and deleted in the loop
    viz::dataprovider::Point thePoint;

	// Loop on the full grid, Y and X first because they are the axis of the plot
	for (PetscInt j = 0; j < My; j++) {
		// Compute y
		y = (double) j * hy;

		for (PetscInt i = 0; i < Mx; i++) {
			// Compute x
			x = (grid[i] + grid[i + 1]) / 2.0 - grid[1];

			// Initialize the value of the concentration to integrate over Z
			double conc = 0.0;

			for (PetscInt k = 0; k < Mz; k++) {
				// If it is the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym && k >= zs
						&& k < zs + zm) {
					// Get the pointer to the beginning of the solution data for this grid point
					gridPointSolution = solutionArray[k][j][i];

					// Integrate over Z
					conc += gridPointSolution[iCluster];
				}
			} // End of the loop on Z

			// Sum all the concentration on Z
			double totalConc = 0.0;
			MPI_Reduce(&conc, &totalConc, 1, MPI_DOUBLE, MPI_SUM, 0,
					PETSC_COMM_WORLD);

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

		// Change the title of the plot and the name of the data
		std::stringstream title;
		title << "First Cluster";
		surfacePlotXY3D->getDataProvider()->setDataName(title.str());
		title << " concentration";
		surfacePlotXY3D->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		surfacePlotXY3D->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		surfacePlotXY3D->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "surfaceXY_TS" << timestep << ".png";
		surfacePlotXY3D->write(fileName.str());
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorSurfaceXZ3D")
/**
 * This is a monitoring method that will save 2D plots of the concentration of
 * a specific cluster at each grid point on the XZ surface, integrating over Y.
 */
PetscErrorCode monitorSurfaceXZ3D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void*) {
	// Initial declarations
	PetscErrorCode ierr;
	const double ****solutionArray, *gridPointSolution;
	PetscInt xs, xm, Mx, ys, ym, My, zs, zm, Mz;
	double x, z;

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
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the physical grid in the x direction
	auto grid = solverHandler.getXGrid();

	// Setup step size variables
	double hz = solverHandler.getStepSizeZ();

	// Choice of the cluster to be plotted
	int iCluster = 0;

	// Create a Point vector to store the data to give to the data provider
	// for the visualization
	auto myPoints = std::make_shared<std::vector<viz::dataprovider::Point> >();
	// Create a point here so that it is not created and deleted in the loop
    viz::dataprovider::Point thePoint;

	// Loop on the full grid, Z and X first because they are the axis of the plot
	for (PetscInt k = 0; k < Mz; k++) {
		// Compute z
		z = (double) k * hz;

		for (PetscInt i = 0; i < Mx; i++) {
			// Compute x
			x = (grid[i] + grid[i + 1]) / 2.0 - grid[1];

			// Initialize the value of the concentration to integrate over Y
			double conc = 0.0;

			for (PetscInt j = 0; j < My; j++) {
				// If it is the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym && k >= zs
						&& k < zs + zm) {
					// Get the pointer to the beginning of the solution data for this grid point
					gridPointSolution = solutionArray[k][j][i];

					// Integrate over Y
					conc += gridPointSolution[iCluster];
				}
			} // End of the loop on Y

			// Sum all the concentration on Y
			double totalConc = 0.0;
			MPI_Reduce(&conc, &totalConc, 1, MPI_DOUBLE, MPI_SUM, 0,
					PETSC_COMM_WORLD);

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

		// Change the title of the plot and the name of the data
		std::stringstream title;
		title << "First Cluster";
		surfacePlotXZ3D->getDataProvider()->setDataName(title.str());
		title << " concentration";
		surfacePlotXZ3D->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		surfacePlotXZ3D->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		surfacePlotXZ3D->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "surfaceXZ_TS" << timestep << ".png";
		surfacePlotXZ3D->write(fileName.str());
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "eventFunction3D")
/**
 * This is a method that checks if the surface should move or bursting happen
 */
PetscErrorCode eventFunction3D(TS ts, PetscReal time, Vec solution,
		PetscScalar *fvalue, void*) {
	// Initial declaration
	PetscErrorCode ierr;
	double ****solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi, Mx, ys, ym, yj, My, zs, zm, zk, Mz;
	fvalue[0] = 1.0, fvalue[1] = 1.0;
	depthPositions3D.clear();

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
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	CHKERRQ(ierr);

	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto &network = solverHandler.getNetwork();

	// Get the physical grid and step size
	auto grid = solverHandler.getXGrid();
	double hy = solverHandler.getStepSizeY();
	double hz = solverHandler.getStepSizeZ();

	// Get the flux handler to know the flux amplitude.
	auto fluxHandler = solverHandler.getFluxHandler();
	double heliumFluxAmplitude = fluxHandler->getFluxAmplitude();

	// Get the delta time from the previous timestep to this timestep
	double dt = time - previousTime;

	// Work of the moving surface first
	if (solverHandler.moveSurface()) {
		// Write the initial surface positions
		if (procId == 0 && util::equal(time, 0.0)) {
			std::ofstream outputFile;
			outputFile.open("surface.txt", std::ios::app);
			outputFile << time << " ";

			// Loop on the possible yj
			for (yj = 0; yj < My; yj++) {
				for (zk = 0; zk < Mz; zk++) {
					// Get the position of the surface at yj, zk
					int surfacePos = solverHandler.getSurfacePosition(yj, zk);
					outputFile << (double) yj * hy << " " << (double) zk * hz
							<< " " << grid[surfacePos + 1] - grid[1] << " ";
				}
			}
			outputFile << std::endl;
			outputFile.close();
		}

		// Get the initial vacancy concentration
		double initialVConc = solverHandler.getInitialVConc();

		// Loop on the possible zk and yj
		for (zk = 0; zk < Mz; zk++) {
			for (yj = 0; yj < My; yj++) {

				// Compute the total density of intersitials that escaped from the
				// surface since last timestep using the stored flux
				nInterstitial3D[yj][zk] += previousIFlux3D[yj][zk] * dt;

				// Remove the sputtering yield since last timestep
				nInterstitial3D[yj][zk] -= sputteringYield3D
						* heliumFluxAmplitude * dt;

				// Get the position of the surface at yj
				int surfacePos = solverHandler.getSurfacePosition(yj, zk);
				xi = surfacePos + solverHandler.getLeftOffset();

				// Initialize the value for the flux
				double newFlux = 0.0;

				// if xi is on this process
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym
						&& zk >= zs && zk < zs + zm) {
					// Get the concentrations at xi = surfacePos + 1
					gridPointSolution = solutionArray[zk][yj][xi];

					// Factor for finite difference
					double hxLeft = 0.0, hxRight = 0.0;
					if (xi - 1 >= 0 && xi < Mx) {
						hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
						hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
					} else if (xi - 1 < 0) {
						hxLeft = (grid[xi + 1] + grid[xi]) / 2.0;
						hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
					} else {
						hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
						hxRight = (grid[xi + 1] - grid[xi]) / 2;
					}
					double factor = 2.0 / (hxLeft + hxRight);

					// Consider each interstitial cluster.
					for (int i = 0; i < iClusterIds3D.size(); i++) {
						auto currId = iClusterIds3D[i];
						// Get the cluster
						auto cluster = network.getClusterCommon(currId);
						// Get its concentration
						double conc = gridPointSolution[currId];
						// Get its size and diffusion coefficient
						int size = i + 1;
						double coef = cluster.getDiffusionCoefficient(xi - xs);
						// Compute the flux going to the left
						newFlux += (double) size * factor * coef * conc
								* hxLeft;
					}
				}

				// Gather newFlux values at this position
				double newTotalFlux = 0.0;
				MPI_Allreduce(&newFlux, &newTotalFlux, 1, MPI_DOUBLE, MPI_SUM,
						PETSC_COMM_WORLD);

				// Update the previous flux
				previousIFlux3D[yj][zk] = newTotalFlux;

				// Compare nInterstitials to the threshold to know if we should move the surface

				// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
				double threshold = (62.8 - initialVConc)
						* (grid[xi] - grid[xi - 1]);
				if (nInterstitial3D[yj][zk] > threshold) {
					// The surface is moving
					fvalue[0] = 0.0;
				}

				// Moving the surface back
				else if (nInterstitial3D[yj][zk] < -threshold / 10.0) {
					// The surface is moving
					fvalue[0] = 0.0;
				}
			}
		}
	}

	// Now work on the bubble bursting
	if (solverHandler.burstBubbles()) {
		using NetworkType =
		core::network::PSIReactionNetwork<core::network::PSIFullSpeciesList>;
		using Spec = typename NetworkType::Species;
		auto psiNetwork = dynamic_cast<NetworkType*>(&network);
		auto dof = network.getDOF();

		// Compute the prefactor for the probability (arbitrary)
		double prefactor = heliumFluxAmplitude * dt * 0.1;

		// The depth parameter to know where the bursting should happen
		double depthParam = solverHandler.getTauBursting();	// nm

		// For now we are not bursting
		bool burst = false;

		// Loop on the full grid
		for (zk = 0; zk < Mz; zk++) {
			for (yj = 0; yj < My; yj++) {
				// Get the surface position
				int surfacePos = solverHandler.getSurfacePosition(yj, zk);
				for (xi = surfacePos + solverHandler.getLeftOffset();
						xi < Mx - solverHandler.getRightOffset(); xi++) {

					// If this is the locally owned part of the grid
					if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym
							&& zk >= zs && zk < zs + zm) {

						// Get the pointer to the beginning of the solution data for this grid point
						gridPointSolution = solutionArray[zk][yj][xi];

						using HostUnmanaged =
						Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
						auto hConcs = HostUnmanaged(gridPointSolution, dof);
						auto dConcs = Kokkos::View<double*>("Concentrations",
								dof);
						deep_copy(dConcs, hConcs);

						// Get the distance from the surface
						double distance = (grid[xi] + grid[xi + 1]) / 2.0
								- grid[surfacePos + 1];

						// Compute the helium density at this grid point
						double heDensity =
								psiNetwork->getTotalAtomConcentration(dConcs,
										Spec::He, 0);

						// Compute the radius of the bubble from the number of helium
						double nV = heDensity * (grid[xi + 1] - grid[xi]) / 4.0;
						//					double nV = pow(heDensity / 5.0, 1.163) * (grid[xi + 1] - grid[xi]);
						double latticeParam = network.getLatticeParameter();
						double tlcCubed = latticeParam * latticeParam
								* latticeParam;
						double radius = (sqrt(3.0) / 4) * latticeParam
								+ cbrt(
										(3.0 * tlcCubed * nV)
												/ (8.0 * core::pi))
								- cbrt(
										(3.0 * tlcCubed)
												/ (8.0 * core::pi));

						// If the radius is larger than the distance to the surface, burst
						if (radius > distance) {
							burst = true;
							depthPositions3D.push_back(
									std::make_tuple(zk, yj, xi));
							// Exit the loop
							continue;
						}
						// Add randomness
						double prob = prefactor
								* (1.0 - (distance - radius) / distance)
								* std::min(1.0,
										exp(
												-(distance - depthParam)
														/ (depthParam * 2.0)));
						double test = solverHandler.getRNG().GetRandomDouble();

						if (prob > test) {
							burst = true;
							depthPositions3D.push_back(
									std::make_tuple(zk, yj, xi));
						}
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
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "postEventFunction3D")
/**
 * This is a method that moves the surface or burst bubbles
 */
PetscErrorCode postEventFunction3D(TS ts, PetscInt nevents,
		PetscInt eventList[], PetscReal time, Vec solution, PetscBool, void*) {

	// Initial declaration
	PetscErrorCode ierr;
	double ****solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi, Mx, ys, ym, yj, My, zs, zm, zk, Mz;

	PetscFunctionBeginUser;

	// Call monitor time hear because it is skipped when post event is used
	ierr = computeFluence(ts, 0, time, solution, NULL);
	CHKERRQ(ierr);
	ierr = monitorTime(ts, 0, time, solution, NULL);
	CHKERRQ(ierr);

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
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	CHKERRQ(ierr);

	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto &network = solverHandler.getNetwork();
	int dof = network.getDOF();

	// Get the physical grid
	auto grid = solverHandler.getXGrid();
	// Get the step sizes
	double hy = solverHandler.getStepSizeY();
	double hz = solverHandler.getStepSizeZ();

	// Take care of bursting
	using NetworkType =
	core::network::PSIReactionNetwork<core::network::PSIFullSpeciesList>;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);

	// Loop on each bursting depth
	for (int i = 0; i < depthPositions3D.size(); i++) {
		// Get the coordinates of the point
		int xi = std::get<2>(depthPositions3D[i]), yj = std::get<1>(
				depthPositions3D[i]), zk = std::get<0>(depthPositions3D[i]);
		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray[zk][yj][xi];

		// Get the surface position
		int surfacePos = solverHandler.getSurfacePosition(yj, zk);
		// Get the distance from the surface
		double distance = (grid[xi] + grid[xi + 1]) / 2.0
				- grid[surfacePos + 1];

		std::cout << "bursting at: " << zk * hz << " " << yj * hy << " "
				<< distance << std::endl;

		// Pinhole case
		// Loop on every cluster
		for (unsigned int i = 0; i < network.getNumClusters(); i++) {
			const auto &clReg =
					psiNetwork->getCluster(i, plsm::onHost).getRegion();
			// Non-grouped clusters
			if (clReg.isSimplex()) {
				// Get the composition
				Composition comp = clReg.getOrigin();
				// Pure He, D, or T case
				if (comp.isOnAxis(Spec::He) || comp.isOnAxis(Spec::D)
						|| comp.isOnAxis(Spec::T)) {
					// Reset concentration
					gridPointSolution[i] = 0.0;
				}
				// Mixed cluster case
				else if (!comp.isOnAxis(Spec::V) && !comp.isOnAxis(Spec::I)) {
					// Transfer concentration to V of the same size
					Composition vComp = Composition::zero();
					vComp[Spec::V] = comp[Spec::V];
					auto vCluster = psiNetwork->findCluster(vComp,
							plsm::onHost);
					gridPointSolution[vCluster.getId()] += gridPointSolution[i];
					gridPointSolution[i] = 0.0;
				}
			}
			// Grouped clusters
			else {
				// Get the factor
				double concFactor = clReg.volume() / clReg[Spec::V].length();
				// Loop on the Vs
				for (auto j : makeIntervalRange(clReg[Spec::V])) {
					// Transfer concentration to V of the same size
					Composition vComp = Composition::zero();
					vComp[Spec::V] = j;
					auto vCluster = psiNetwork->findCluster(vComp,
							plsm::onHost);
					// TODO: refine formula with V moment
					gridPointSolution[vCluster.getId()] += gridPointSolution[i]
							* concFactor;
				}

				// Reset the concentration and moments
				gridPointSolution[i] = 0.0;
				auto momentIds =
						psiNetwork->getCluster(i, plsm::onHost).getMomentIds();
				for (std::size_t j = 0; j < momentIds.extent(0); j++) {
					gridPointSolution[momentIds(j)] = 0.0;
				}
			}
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

	// Loop on the possible zk and yj
	for (zk = 0; zk < Mz; zk++) {
		for (yj = 0; yj < My; yj++) {
			// Get the position of the surface at yj
			int surfacePos = solverHandler.getSurfacePosition(yj, zk);
			xi = surfacePos + solverHandler.getLeftOffset();

			// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
			double threshold = (62.8 - initialVConc)
					* (grid[xi] - grid[xi - 1]);

			// Move the surface up
			if (nInterstitial3D[yj][zk] > threshold) {
				int nGridPoints = 0;
				// Move the surface up until it is smaller than the next threshold
				while (nInterstitial3D[yj][zk] > threshold) {
					// Move the surface higher
					surfacePos--;
					xi = surfacePos + solverHandler.getLeftOffset();
					nGridPoints++;
					// Update the number of interstitials
					nInterstitial3D[yj][zk] -= threshold;
					// Update the thresold
					double threshold = (62.8 - initialVConc)
							* (grid[xi] - grid[xi - 1]);
				}

				// Throw an exception if the position is negative
				if (surfacePos < 0) {
					PetscBool flagCheck;
					ierr = PetscOptionsHasName(NULL, NULL, "-check_collapse",
							&flagCheck);
					CHKERRQ(ierr);
					if (flagCheck) {
						// Write the convergence reason
						std::ofstream outputFile;
						outputFile.open("solverStatus.txt");
						outputFile << "overgrid" << std::endl;
						outputFile.close();
					}
					throw std::string(
							"\nxolotlSolver::Monitor3D: The surface is trying to go outside of the grid!!");
				}

				// Printing information about the extension of the material
				if (procId == 0) {
					std::cout << "Adding " << nGridPoints
							<< " points to the grid on " << yj * hy << " and "
							<< zk * hz << " at time: " << time << " s."
							<< std::endl;
				}

				// Set it in the solver
				solverHandler.setSurfacePosition(surfacePos, yj, zk);

				// Initialize the vacancy concentration and the temperature on the new grid points
				// Get the single vacancy ID
				auto singleVacancyCluster = network.getSingleVacancy();
				auto vacancyIndex =
						core::network::IReactionNetwork::invalidIndex();
				if (singleVacancyCluster.getId()
						!= core::network::IReactionNetwork::invalidIndex())
					vacancyIndex = singleVacancyCluster.getId();
				// Get the surface temperature
				double temp = 0.0;
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym
						&& zk >= zs && zk < zs + zm) {
					temp = solutionArray[zk][yj][xi][dof];
				}
				double surfTemp = 0.0;
				MPI_Allreduce(&temp, &surfTemp, 1, MPI_DOUBLE, MPI_SUM,
						PETSC_COMM_WORLD);
				// Loop on the new grid points
				while (nGridPoints >= 0) {
					// Position of the newly created grid point
					xi = surfacePos + nGridPoints;

					// If xi is on this process
					if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym
							&& zk >= zs && zk < zs + zm) {
						// Get the concentrations
						gridPointSolution = solutionArray[zk][yj][xi];

						// Set the new surface temperature
						gridPointSolution[dof] = surfTemp;

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
			else if (nInterstitial3D[yj][zk] < -threshold / 10.0) {
				// Move it back as long as the number of interstitials in negative
				while (nInterstitial3D[yj][zk] < 0.0) {
					// Compute the threshold to a deeper grid point
					threshold = (62.8 - initialVConc)
							* (grid[xi + 1] - grid[xi]);
					// Set all the concentrations to 0.0 at xi = surfacePos + 1
					// if xi is on this process
					if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym
							&& zk >= zs && zk < zs + zm) {
						// Get the concentrations at xi = surfacePos + 1
						gridPointSolution = solutionArray[zk][yj][xi];
						// Loop on DOF
						for (int i = 0; i < dof; i++) {
							gridPointSolution[i] = 0.0;
						}
					}

					// Move the surface deeper
					surfacePos++;
					xi = surfacePos + solverHandler.getLeftOffset();
					// Update the number of interstitials
					nInterstitial3D[yj][zk] += threshold;
				}

				// Printing information about the extension of the material
				if (procId == 0) {
					std::cout << "Removing grid points to the grid on "
							<< yj * hy << " and " << zk * hz << " at time: "
							<< time << " s." << std::endl;
				}

				// Set it in the solver
				solverHandler.setSurfacePosition(surfacePos, yj, zk);
			}
		}
	}
	// Get the modified trap-mutation handler to reinitialize it
	auto mutationHandler = solverHandler.getMutationHandler();
	auto advecHandlers = solverHandler.getAdvectionHandlers();

	// Get the vector of positions of the surface
	std::vector<std::vector<int> > surfaceIndices;
	for (PetscInt i = 0; i < My; i++) {
		// Create a temporary vector
		std::vector<int> temp;
		for (PetscInt j = 0; j < Mz; j++) {
			temp.push_back(solverHandler.getSurfacePosition(i, j));
		}
		// Add the temporary vector to the vector of surface indices
		surfaceIndices.push_back(temp);
	}

	mutationHandler->initializeIndex3D(surfaceIndices,
			solverHandler.getNetwork(), advecHandlers, grid, xm, xs, ym, hy, ys,
			zm, hz, zs);

	// Write the surface positions
	if (procId == 0) {
		std::ofstream outputFile;
		outputFile.open("surface.txt", std::ios::app);
		outputFile << time << " ";

		// Loop on the possible yj
		for (yj = 0; yj < My; yj++) {
			for (zk = 0; zk < Mz; zk++) {
				// Get the position of the surface at yj, zk
				int surfacePos = solverHandler.getSurfacePosition(yj, zk);
				outputFile << (double) yj * hy << " " << (double) zk * hz << " "
						<< grid[surfacePos + 1] - grid[1] << " ";
			}
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
PetscErrorCode setupPetsc3DMonitor(TS ts) {
	PetscErrorCode ierr;

	// Get the process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the xolotlViz handler registry
	auto vizHandlerRegistry = factory::viz::getVizHandlerRegistry();

	// Flags to launch the monitors or not
	PetscBool flagCheck, flagPerf, flagHeRetention, flagXeRetention, flagStatus,
			flag2DXYPlot, flag2DXZPlot;

	// Check the option -check_collapse
	ierr = PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCheck);
	checkPetscError(ierr,
			"setupPetsc3DMonitor: PetscOptionsHasName (-check_collapse) failed.");

	// Check the option -plot_perf
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_perf", &flagPerf);
	checkPetscError(ierr,
			"setupPetsc3DMonitor: PetscOptionsHasName (-plot_perf) failed.");

	// Check the option -plot_2d_xy
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_2d_xy", &flag2DXYPlot);
	checkPetscError(ierr,
			"setupPetsc3DMonitor: PetscOptionsHasName (-plot_2d_xy) failed.");

	// Check the option -plot_2d_xz
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_2d_xz", &flag2DXZPlot);
	checkPetscError(ierr,
			"setupPetsc3DMonitor: PetscOptionsHasName (-plot_2d_xz) failed.");

	// Check the option -helium_retention
	ierr = PetscOptionsHasName(NULL, NULL, "-helium_retention",
			&flagHeRetention);
	checkPetscError(ierr,
			"setupPetsc3DMonitor: PetscOptionsHasName (-helium_retention) failed.");

	// Check the option -xenon_retention
	ierr = PetscOptionsHasName(NULL, NULL, "-xenon_retention",
			&flagXeRetention);
	checkPetscError(ierr,
			"setupPetsc3DMonitor: PetscOptionsHasName (-xenon_retention) failed.");

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, NULL, "-start_stop", &flagStatus);
	checkPetscError(ierr,
			"setupPetsc3DMonitor: PetscOptionsHasName (-start_stop) failed.");

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto &network = solverHandler.getNetwork();

	// Determine if we have an existing restart file,
	// and if so, it it has had timesteps written to it.
	std::unique_ptr<io::XFile> networkFile;
	std::unique_ptr<io::XFile::TimestepGroup> lastTsGroup;
	std::string networkName = solverHandler.getNetworkName();
	bool hasConcentrations = false;
	if (not networkName.empty()) {
		networkFile.reset(new io::XFile(networkName));
		auto concGroup = networkFile->getGroup<
				io::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
		if (hasConcentrations) {
			lastTsGroup = concGroup->getLastTimestepGroup();
		}
	}

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);
	checkPetscError(ierr, "setupPetsc3DMonitor: TSGetDM failed.");

	// Get the total size of the grid
	PetscInt Mx, My, Mz;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);
	checkPetscError(ierr, "setupPetsc3DMonitor: DMDAGetInfo failed.");

	// Set the post step processing to stop the solver if the time step collapses
	if (flagCheck) {
		// Find the threshold
		PetscBool flag;
		ierr = PetscOptionsGetReal(NULL, NULL, "-check_collapse",
				&timeStepThreshold, &flag);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: PetscOptionsGetInt (-check_collapse) failed.");
		if (!flag)
			timeStepThreshold = 1.0e-16;

		// Set the post step process that tells the solver when to stop if the time step collapse
		ierr = TSSetPostStep(ts, checkTimeStep);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: TSSetPostStep (checkTimeStep) failed.");
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr = PetscOptionsGetReal(NULL, NULL, "-start_stop", &hdf5Stride3D,
				&flag);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: PetscOptionsGetInt (-start_stop) failed.");
		if (!flag)
			hdf5Stride3D = 1.0;

		// Compute the correct hdf5Previous3D for a restart
		if (hasConcentrations) {
			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			previousTime = lastTsGroup->readPreviousTime();
			hdf5Previous3D = (int) (previousTime / hdf5Stride3D);
		}

		// Don't do anything if both files have the same name
		if (hdf5OutputName3D != solverHandler.getNetworkName()) {

			// Get the solver handler
			auto &solverHandler = PetscSolver::getSolverHandler();

			// Get the physical grid in the x direction
			auto grid = solverHandler.getXGrid();

			// Setup step size variables
			double hy = solverHandler.getStepSizeY();
			double hz = solverHandler.getStepSizeZ();

			// Create a checkpoint file.
			// Create and initialize a checkpoint file.
			// We do this in its own scope so that the file
			// is closed when the file object goes out of scope.
			// We want it to close before we (potentially) copy
			// the network from another file using a single-process
			// MPI communicator.
			{
				io::XFile checkpointFile(hdf5OutputName3D, grid,
						PETSC_COMM_WORLD, My, hy, Mz, hz);
			}

			// Copy the network group from the given file (if it has one).
			// We open the files using a single-process MPI communicator
			// because it is faster for a single process to do the
			// copy with HDF5's H5Ocopy implementation than it is
			// when all processes call the copy function.
			// The checkpoint file must be closed before doing this.
			writeNetwork(PETSC_COMM_WORLD, solverHandler.getNetworkName(),
					hdf5OutputName3D, network);
		}

		// startStop3D will be called at each timestep
		ierr = TSMonitorSet(ts, startStop3D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: TSMonitorSet (startStop3D) failed.");
	}

	// If the user wants the surface to be able to move
	if (solverHandler.moveSurface() || solverHandler.burstBubbles()) {
		// Surface
		if (solverHandler.moveSurface()) {
			// Initialize nInterstitial3D and previousIFlux3D before monitoring the
			// interstitial flux
			for (PetscInt j = 0; j < My; j++) {
				// Create a one dimensional vector of double
				std::vector<double> tempVector;
				for (PetscInt k = 0; k < Mz; k++) {
					tempVector.push_back(0.0);
				}
				// Add the tempVector to nInterstitial3D and previousIFlux3D
				// to create their initial structure
				nInterstitial3D.push_back(tempVector);
				previousIFlux3D.push_back(tempVector);
			}

			// Get the interstitial information at the surface if concentrations were stored
			if (hasConcentrations) {
				// Get the interstitial quantity from the HDF5 file
				nInterstitial3D = lastTsGroup->readData3D("nInterstitial");
				// Get the previous I flux from the HDF5 file
				previousIFlux3D = lastTsGroup->readData3D("previousIFlux");
				// Get the previous time from the HDF5 file
				previousTime = lastTsGroup->readPreviousTime();
			}

			// Get the sputtering yield
			sputteringYield3D = solverHandler.getSputteringYield();

			// Master process
			if (procId == 0) {
				// Clear the file where the surface will be written
				std::ofstream outputFile;
				outputFile.open("surface.txt");
				outputFile.close();
			}
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
		ierr = TSSetEventHandler(ts, 2, direction, terminate, eventFunction3D,
				postEventFunction3D, NULL);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: TSSetEventHandler (eventFunction3D) failed.");
	}

	// Set the monitor to save performance plots (has to be in parallel)
	if (flagPerf) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a ScatterPlot
			perfPlot = vizHandlerRegistry->getPlot("perfPlot",
					viz::PlotType::SCATTER);

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>(
					"labelProvider");
			labelProvider->axis1Label = "Process ID";
			labelProvider->axis2Label = "Solver Time";

			// Give it to the plot
			perfPlot->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider = std::make_shared<viz::dataprovider::CvsXDataProvider>(
					"dataProvider");

			// Give it to the plot
			perfPlot->setDataProvider(dataProvider);
		}

		// monitorPerf will be called at each timestep
		ierr = TSMonitorSet(ts, monitorPerf, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: TSMonitorSet (monitorPerf) failed.");
	}

	// Set the monitor to compute the helium fluence for the retention calculation
	if (flagHeRetention) {

		// Get the previous time if concentrations were stored and initialize the fluence
		if (hasConcentrations) {
			// Get the previous time from the HDF5 file
			double time = lastTsGroup->readPreviousTime();
			// Initialize the fluence
			auto fluxHandler = solverHandler.getFluxHandler();
			// The length of the time step
			double dt = time;
			// Increment the fluence with the value at this current timestep
			fluxHandler->incrementFluence(dt);
			// Get the previous time from the HDF5 file
			// TODO is this same as 'time' above?
			previousTime = lastTsGroup->readPreviousTime();
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeFluence, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeHeliumRetention3D will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumRetention3D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: TSMonitorSet (computeHeliumRetention3D) failed.");

		// Master process
		if (procId == 0) {
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("retentionOut.txt");
			outputFile.close();
		}
	}

	// Set the monitor to compute the xenon fluence and the retention
	// for the retention calculation
	if (flagXeRetention) {
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
			// TODO isn't this the same as 'time' above?
			previousTime = lastTsGroup->readPreviousTime();
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeFluence, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeXenonRetention3D will be called at each timestep
		ierr = TSMonitorSet(ts, computeXenonRetention3D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: TSMonitorSet (computeXenonRetention3D) failed.");

		// Master process
		if (procId == 0) {
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("retentionOut.txt");
			outputFile.close();
		}
	}

	// Set the monitor to save surface plots of clusters concentration
	if (flag2DXYPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a SurfacePlot
			surfacePlotXY3D = vizHandlerRegistry->getPlot("surfacePlotXY3D",
					viz::PlotType::SURFACE);

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>(
					"labelProvider");
			labelProvider->axis1Label = "Depth (nm)";
			labelProvider->axis2Label = "Y (nm)";
			labelProvider->axis3Label = "Concentration";

			// Give it to the plot
			surfacePlotXY3D->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider = std::make_shared<viz::dataprovider::CvsXYDataProvider>(
					"dataProvider");

			// Give it to the plot
			surfacePlotXY3D->setDataProvider(dataProvider);
		}

		// monitorSurfaceXY3D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSurfaceXY3D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: TSMonitorSet (monitorSurfaceXY3D) failed.");
	}

	// Set the monitor to save surface plots of clusters concentration
	if (flag2DXZPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a SurfacePlot
			surfacePlotXZ3D = vizHandlerRegistry->getPlot("surfacePlotXZ3D",
					viz::PlotType::SURFACE);

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>(
					"labelProvider");
			labelProvider->axis1Label = "Depth (nm)";
			labelProvider->axis2Label = "Z (nm)";
			labelProvider->axis3Label = "Concentration";

			// Give it to the plot
			surfacePlotXZ3D->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider = std::make_shared<viz::dataprovider::CvsXYDataProvider>(
					"dataProvider");

			// Give it to the plot
			surfacePlotXZ3D->setDataProvider(dataProvider);
		}

		// monitorSurfaceXZ3D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSurfaceXZ3D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc3DMonitor: TSMonitorSet (monitorSurfaceXZ3D) failed.");
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(ts, monitorTime, NULL, NULL);
	checkPetscError(ierr,
			"setupPetsc3DMonitor: TSMonitorSet (monitorTime) failed.");

	PetscFunctionReturn(0);
}

} /* end namespace monitor */
} /* end namespace solver */
} /* end namespace xolotl */
