// Includes
#include <petscsys.h>
#include <petscts.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/io/XFile.h>
#include <xolotl/solver/PetscSolver.h>
#include <xolotl/solver/monitor/Monitor.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/MathUtils.h>
#include <xolotl/util/RandomNumberGenerator.h>
#include <xolotl/viz/IPlot.h>
#include <xolotl/viz/LabelProvider.h>
#include <xolotl/viz/PlotType.h>
#include <xolotl/viz/VizHandlerRegistry.h>
#include <xolotl/viz/dataprovider/CvsXDataProvider.h>
#include <xolotl/viz/dataprovider/CvsXYDataProvider.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
// Declaration of the functions defined in Monitor.cpp
extern PetscErrorCode
checkTimeStep(TS ts);
extern PetscErrorCode
monitorTime(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);
extern PetscErrorCode
computeFluence(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);
extern PetscErrorCode
monitorPerf(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);

// Declaration of the variables defined in Monitor.cpp
extern std::shared_ptr<viz::IPlot> perfPlot;
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
std::vector<std::vector<double>> previousIFlux3D;
//! The variable to store the total number of interstitials going through the
//! surface.
std::vector<std::vector<double>> nInterstitial3D;
//! The variable to store the sputtering yield at the surface.
double sputteringYield3D = 0.0;
// The vector of depths at which bursting happens
std::vector<std::array<PetscInt, 3>> depthPositions3D;
// The vector of ids for diffusing interstitial clusters
std::vector<IdType> iClusterIds3D;
// The id of the largest cluster
int largestClusterId3D = -1;
// The concentration threshold for the largest cluster
double largestThreshold3D = 1.0e-12;
// Tracks the previous TS number
PetscInt previousTSNumber3D = -1;

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorLargest3D")
/**
 * This is a monitoring method that looks at the largest cluster concentration
 */
PetscErrorCode
monitorLargest3D(TS ts, PetscInt timestep, PetscReal time, Vec solution, void*)
{
	// Initial declaration
	PetscErrorCode ierr;
	double ****solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the solver handler and local coordinates
	auto& solverHandler = PetscSolver::getSolverHandler();
	solverHandler.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Loop on the local grid
	for (auto k = zs; k < zs + zm; k++)
		for (auto j = ys; j < ys + ym; j++)
			for (auto i = xs; i < xs + xm; i++) {
				// Get the pointer to the beginning of the solution data for
				// this grid point
				gridPointSolution = solutionArray[k][j][i];
				// Check the concentration
				if (gridPointSolution[largestClusterId3D] >
					largestThreshold3D) {
					ierr = TSSetConvergedReason(ts, TS_CONVERGED_USER);
					CHKERRQ(ierr);
					// Send an error
					throw std::runtime_error(
						"\nxolotlSolver::Monitor3D: The largest "
						"cluster concentration is too high!!");
				}
			}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "startStop3D")
/**
 * This is a monitoring method that will update an hdf5 file every given time.
 */
PetscErrorCode
startStop3D(TS ts, PetscInt timestep, PetscReal time, Vec solution, void*)
{
	// Initial declarations
	PetscErrorCode ierr;
	const double ****solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Get the solver handler and local coordinates
	auto& solverHandler = PetscSolver::getSolverHandler();
	solverHandler.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Compute the dt
	double previousTime = solverHandler.getPreviousTime();
	double dt = time - previousTime;

	// Don't do anything if it is not on the stride
	if (((PetscInt)((time + dt / 10.0) / hdf5Stride3D) <= hdf5Previous3D) &&
		timestep > 0)
		PetscFunctionReturn(0);

	// Update the previous time
	if ((PetscInt)((time + dt / 10.0) / hdf5Stride3D) > hdf5Previous3D)
		hdf5Previous3D++;

	// Gets the process ID (important when it is running in parallel)
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the network and dof
	auto& network = solverHandler.getNetwork();
	const auto dof = network.getDOF();

	// Create an array for the concentration
	double concArray[dof + 1][2];

	// Get the vector of positions of the surface
	std::vector<std::vector<int>> surfaceIndices;
	for (auto i = 0; i < My; i++) {
		// Create a temporary vector
		std::vector<int> temp;
		for (auto j = 0; j < Mz; j++) {
			temp.push_back(solverHandler.getSurfacePosition(i, j));
		}
		// Add the temporary vector to the vector of surface indices
		surfaceIndices.push_back(temp);
	}

	// Open the existing HDF5 file.
	io::XFile checkpointFile(
		hdf5OutputName3D, xolotlComm, io::XFile::AccessMode::OpenReadWrite);

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	CHKERRQ(ierr);

	// Add a concentration sub group
	auto concGroup = checkpointFile.getGroup<io::XFile::ConcentrationGroup>();
	assert(concGroup);
	auto tsGroup = concGroup->addTimestepGroup(
		timestep, time, previousTime, currentTimeStep);

	if (solverHandler.moveSurface()) {
		// Write the surface positions in the concentration sub group
		tsGroup->writeSurface3D(
			surfaceIndices, nInterstitial3D, previousIFlux3D);
	}

	// Loop on the full grid
	for (auto k = 0; k < Mz; k++) {
		for (auto j = 0; j < My; j++) {
			for (auto i = 0; i < Mx; i++) {
				// Wait for all the processes
				MPI_Barrier(xolotlComm);

				// Size of the concentration that will be stored
				int concSize = -1;
				// To save which proc has the information
				int concId = 0;
				// To know which process should write
				bool write = false;

				// If it is the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym &&
					k >= zs && k < zs + zm) {
					write = true;
					// Get the pointer to the beginning of the solution data for
					// this grid point
					gridPointSolution = solutionArray[k][j][i];

					// Loop on the concentrations
					for (auto l = 0; l < dof + 1; l++) {
						if (std::fabs(gridPointSolution[l]) > 1.0e-16) {
							// Increase concSize
							concSize++;
							// Fill the concArray
							concArray[concSize][0] = (double)l;
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
				MPI_Allreduce(
					&concId, &concProc, 1, MPI_INT, MPI_SUM, xolotlComm);

				// Broadcast the size
				MPI_Bcast(&concSize, 1, MPI_INT, concProc, xolotlComm);

				// Skip the grid point if the size is 0
				if (concSize == 0)
					continue;

				// All processes create the dataset and fill it
				tsGroup->writeConcentrationDataset(
					concSize, concArray, write, i, j, k);
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
PetscErrorCode
computeHeliumRetention3D(TS ts, PetscInt, PetscReal time, Vec solution, void*)
{
	// Initial declarations
	PetscErrorCode ierr;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Get the solver handler and local coordinates
	auto& solverHandler = PetscSolver::getSolverHandler();
	solverHandler.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = solverHandler.getFluxHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the physical grid in the x direction
	auto grid = solverHandler.getXGrid();

	// Get the network
	using NetworkType = core::network::IPSIReactionNetwork;
	auto& network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	const auto dof = network.getDOF();

	// Setup step size variables
	double hy = solverHandler.getStepSizeY();
	double hz = solverHandler.getStepSizeZ();

	// Get the array of concentration
	double ****solutionArray, *gridPointSolution;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration over the grid
	auto numSpecies = network.getSpeciesListSize();
	auto myConcData = std::vector<double>(numSpecies, 0.0);

	// Loop on the grid
	for (auto zk = zs; zk < zs + zm; zk++) {
		for (auto yj = ys; yj < ys + ym; yj++) {
			// Get the surface position
			auto surfacePos = solverHandler.getSurfacePosition(yj, zk);

			for (auto xi = xs; xi < xs + xm; xi++) {
				// Boundary conditions
				if (xi < surfacePos + solverHandler.getLeftOffset() ||
					xi >= Mx - solverHandler.getRightOffset())
					continue;

				// Get the pointer to the beginning of the solution data for
				// this grid point
				gridPointSolution = solutionArray[zk][yj][xi];

				double hx = grid[xi + 1] - grid[xi];

				using HostUnmanaged = Kokkos::View<double*, Kokkos::HostSpace,
					Kokkos::MemoryUnmanaged>;
				auto hConcs = HostUnmanaged(gridPointSolution, dof);
				auto dConcs = Kokkos::View<double*>("Concentrations", dof);
				deep_copy(dConcs, hConcs);

				// Get the total concentrations at this grid point
				for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
					myConcData[id()] +=
						network.getTotalAtomConcentration(dConcs, id, 1) * hx *
						hy * hz;
				}
			}
		}
	}

	// Get the current process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Sum all the concentrations through MPI reduce
	auto totalConcData = std::vector<double>(numSpecies, 0.0);

	MPI_Reduce(myConcData.data(), totalConcData.data(), myConcData.size(),
		MPI_DOUBLE, MPI_SUM, 0, xolotlComm);

	// Master process
	if (procId == 0) {
		// Get the total size of the grid rescale the concentrations
		PetscInt Mx, My, Mz;
		ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
		CHKERRQ(ierr);

		// Compute the total surface irradiated by the helium flux
		double surface = (double)(My * Mz) * hy * hz;

		// Rescale the concentration
		for (auto i = 0; i < numSpecies; ++i) {
			totalConcData[i] /= surface;
		}

		// Get the fluence
		double fluence = fluxHandler->getFluence();

		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			std::cout << network.getSpeciesName(id)
					  << " content = " << totalConcData[id()] << '\n';
		}
		std::cout << "Fluence = " << fluence << "\n" << std::endl;

		// Uncomment to write the retention and the fluence in a file
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt", std::ios::app);
		outputFile << fluence << " ";
		for (auto i = 0; i < numSpecies; ++i) {
			outputFile << totalConcData[i] << " ";
		}
		outputFile << std::endl;
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
PetscErrorCode
computeXenonRetention3D(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void*)
{
	// Initial declarations
	PetscErrorCode ierr;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Get the solver handler and local coordinates
	auto& solverHandler = PetscSolver::getSolverHandler();
	solverHandler.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// Setup step size variables
	double hy = solverHandler.getStepSizeY();
	double hz = solverHandler.getStepSizeZ();

	using NetworkType = core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Degrees of freedom is the total number of clusters in the network
	auto& network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	const auto dof = network.getDOF();

	// Get the complete data array, including ghost cells
	Vec localSolution;
	ierr = DMGetLocalVector(da, &localSolution);
	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da, solution, INSERT_VALUES, localSolution);
	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, solution, INSERT_VALUES, localSolution);
	CHKERRQ(ierr);
	// Get the array of concentration
	PetscReal ****solutionArray, *gridPointSolution;
	ierr = DMDAVecGetArrayDOFRead(da, localSolution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration and other values over the grid
	double xeConcentration = 0.0, bubbleConcentration = 0.0, radii = 0.0,
		   partialBubbleConcentration = 0.0, partialRadii = 0.0;

	// Get the minimum size for the radius
	auto minSizes = solverHandler.getMinSizes();

	// Get Xe_1
	Composition xeComp = Composition::zero();
	xeComp[Spec::Xe] = 1;
	auto xeCluster = network.findCluster(xeComp, plsm::onHost);
	auto xeId = xeCluster.getId();

	// Loop on the grid
	for (auto zk = zs; zk < zs + zm; zk++) {
		for (auto yj = ys; yj < ys + ym; yj++) {
			for (auto xi = xs; xi < xs + xm; xi++) {
				// Get the pointer to the beginning of the solution data for
				// this grid point
				gridPointSolution = solutionArray[zk][yj][xi];

				using HostUnmanaged = Kokkos::View<double*, Kokkos::HostSpace,
					Kokkos::MemoryUnmanaged>;
				auto hConcs = HostUnmanaged(gridPointSolution, dof);
				auto dConcs = Kokkos::View<double*>("Concentrations", dof);
				deep_copy(dConcs, hConcs);

				double hx = grid[xi + 1] - grid[xi];

				// Get the concentrations
				xeConcentration +=
					network.getTotalAtomConcentration(dConcs, Spec::Xe, 1) *
					hx * hy * hz;
				bubbleConcentration +=
					network.getTotalConcentration(dConcs, Spec::Xe, 1) * hx *
					hy * hz;
				radii +=
					network.getTotalRadiusConcentration(dConcs, Spec::Xe, 1) *
					hx * hy * hz;
				partialBubbleConcentration +=
					network.getTotalConcentration(
						dConcs, Spec::Xe, minSizes[0]) *
					hx * hy * hz;
				partialRadii += network.getTotalRadiusConcentration(
									dConcs, Spec::Xe, minSizes[0]) *
					hx * hy * hz;

				// Set the volume fraction
				double volumeFrac = network.getTotalVolumeFraction(
					dConcs, Spec::Xe, minSizes[0]);
				solverHandler.setVolumeFraction(
					volumeFrac, xi - xs, yj - ys, zk - zs);
				solverHandler.setMonomerConc(
					gridPointSolution[xeCluster.getId()], xi - xs, yj - ys,
					zk - zs);
			}
		}
	}

	// Get the current process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Sum all the concentrations through MPI reduce
	std::array<double, 5> myConcData{xeConcentration, bubbleConcentration,
		radii, partialBubbleConcentration, partialRadii};
	std::array<double, 5> totalConcData{0.0, 0.0, 0.0, 0.0, 0.0};
	MPI_Reduce(myConcData.data(), totalConcData.data(), myConcData.size(),
		MPI_DOUBLE, MPI_SUM, 0, xolotlComm);

	// GB
	// Get the delta time from the previous timestep to this timestep
	double dt = time - solverHandler.getPreviousTime();
	// Sum and gather the previous flux
	double globalXeFlux = 0.0;
	// Get the vector from the solver handler
	auto gbVector = solverHandler.getGBVector();
	// Get the previous flux vector
	auto& localNE = solverHandler.getLocalNE();
	// Loop on the GB
	for (auto const& pair : gbVector) {
		// Middle
		auto xi = std::get<0>(pair);
		auto yj = std::get<1>(pair);
		auto zk = std::get<2>(pair);
		// Check we are on the right proc
		if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym && zk >= zs &&
			zk < zs + zm) {
			double previousXeFlux =
				std::get<1>(localNE[xi - xs][yj - ys][zk - zs]);
			globalXeFlux +=
				previousXeFlux * (grid[xi + 1] - grid[xi]) * hy * hz;
			// Set the amount in the vector we keep
			solverHandler.setLocalXeRate(
				previousXeFlux * dt, xi - xs, yj - ys, zk - zs);
		}
	}
	double totalXeFlux = 0.0;
	MPI_Reduce(
		&globalXeFlux, &totalXeFlux, 1, MPI_DOUBLE, MPI_SUM, 0, xolotlComm);
	// Master process
	if (procId == 0) {
		// Get the previous value of Xe that went to the GB
		double nXenon = solverHandler.getNXeGB();
		// Compute the total number of Xe that went to the GB
		nXenon += totalXeFlux * dt;
		solverHandler.setNXeGB(nXenon);
	}

	// Get the number of species
	auto numSpecies = network.getSpeciesListSize();

	// Get the vector of diffusing clusters
	auto diffusionHandler = solverHandler.getDiffusionHandler();
	auto diffusingIds = diffusionHandler->getDiffusingIds();

	// Loop on the GB
	for (auto const& pair : gbVector) {
		// Local rate
		auto myRate = std::vector<double>(numSpecies, 0.0);
		// Define left and right with reference to the middle point
		// Middle
		auto xi = std::get<0>(pair);
		auto yj = std::get<1>(pair);
		auto zk = std::get<2>(pair);

		// Factor for finite difference
		double hxLeft = 0.0, hxRight = 0.0;
		if (xi >= 1 && xi < Mx) {
			hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
			hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
		}
		else if (xi < 1) {
			hxLeft = grid[xi + 1] - grid[xi];
			hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
		}
		else {
			hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
			hxRight = grid[xi + 1] - grid[xi];
		}
		double factor = 2.0 / (hxLeft + hxRight);

		// Check we are on the right proc
		if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym && zk >= zs &&
			zk < zs + zm) {
			// X segment
			// Left
			xi = std::get<0>(pair) - 1;
			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[zk][yj][xi];
			// Compute the flux coming from the left
			network.updateOutgoingDiffFluxes(gridPointSolution, factor / hxLeft,
				diffusingIds, myRate, xi + 1 - xs);

			// Right
			xi = std::get<0>(pair) + 1;
			gridPointSolution = solutionArray[zk][yj][xi];
			// Compute the flux coming from the right
			network.updateOutgoingDiffFluxes(gridPointSolution,
				factor / hxRight, diffusingIds, myRate, xi + 1 - xs);

			// Y segment
			// Bottom
			xi = std::get<0>(pair);
			yj = std::get<1>(pair) - 1;
			gridPointSolution = solutionArray[zk][yj][xi];
			// Compute the flux coming from the bottom
			network.updateOutgoingDiffFluxes(gridPointSolution, 1.0 / (hy * hy),
				diffusingIds, myRate, xi + 1 - xs);

			// Top
			yj = std::get<1>(pair) + 1;
			gridPointSolution = solutionArray[zk][yj][xi];
			// Compute the flux coming from the top
			network.updateOutgoingDiffFluxes(gridPointSolution, 1.0 / (hy * hy),
				diffusingIds, myRate, xi + 1 - xs);

			// Z segment
			// Back
			yj = std::get<1>(pair);
			zk = std::get<2>(pair) - 1;
			gridPointSolution = solutionArray[zk][yj][xi];
			// Compute the flux coming from the back
			network.updateOutgoingDiffFluxes(gridPointSolution, 1.0 / (hz * hz),
				diffusingIds, myRate, xi + 1 - xs);

			// Front
			zk = std::get<2>(pair) + 1;
			gridPointSolution = solutionArray[zk][yj][xi];
			// Compute the flux coming from the front
			network.updateOutgoingDiffFluxes(gridPointSolution, 1.0 / (hz * hz),
				diffusingIds, myRate, xi + 1 - xs);

			// Middle
			xi = std::get<0>(pair);
			yj = std::get<1>(pair);
			zk = std::get<2>(pair);
			solverHandler.setPreviousXeFlux(
				myRate[0], xi - xs, yj - ys, zk - zs);
		}
	}

	// Master process
	if (procId == 0) {
		// Compute the total surface irradiated
		double surface = (double)My * hy * (double)Mz * hz;
		// Get the number of Xe that went to the GB
		double nXenon = solverHandler.getNXeGB();

		totalConcData[0] = totalConcData[0] / surface;

		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		std::cout << "Xenon concentration = " << totalConcData[0] << std::endl;
		std::cout << "Xenon GB = " << nXenon / surface << std::endl
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
				   << averagePartialRadius << " " << nXenon / surface
				   << std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, localSolution, &solutionArray);
	CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da, &localSolution);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorSurfaceXY3D")
/**
 * This is a monitoring method that will save 2D plots of the concentration of
 * a specific cluster at each grid point on the XY surface, integrating over Z.
 */
PetscErrorCode
monitorSurfaceXY3D(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void*)
{
	// Initial declarations
	PetscErrorCode ierr;
	const double ****solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;
	double x, y;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % 10 != 0)
		PetscFunctionReturn(0);

	// Gets the process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the solver handler and local coordinates
	auto& solverHandler = PetscSolver::getSolverHandler();
	solverHandler.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the physical grid in the x direction
	auto grid = solverHandler.getXGrid();

	// Setup step size variables
	double hy = solverHandler.getStepSizeY();

	// Choice of the cluster to be plotted
	IdType iCluster = 0;

	// Create a DataPoint vector to store the data to give to the data provider
	// for the visualization
	auto myPoints =
		std::make_shared<std::vector<viz::dataprovider::DataPoint>>();
	// Create a point here so that it is not created and deleted in the loop
	viz::dataprovider::DataPoint thePoint;

	// Loop on the full grid, Y and X first because they are the axis of the
	// plot
	for (auto j = 0; j < My; j++) {
		// Compute y
		y = (double)j * hy;

		for (auto i = 0; i < Mx; i++) {
			// Compute x
			x = (grid[i] + grid[i + 1]) / 2.0 - grid[1];

			// Initialize the value of the concentration to integrate over Z
			double conc = 0.0;

			for (auto k = 0; k < Mz; k++) {
				// If it is the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym &&
					k >= zs && k < zs + zm) {
					// Get the pointer to the beginning of the solution data for
					// this grid point
					gridPointSolution = solutionArray[k][j][i];

					// Integrate over Z
					conc += gridPointSolution[iCluster];
				}
			} // End of the loop on Z

			// Sum all the concentration on Z
			double totalConc = 0.0;
			MPI_Reduce(
				&conc, &totalConc, 1, MPI_DOUBLE, MPI_SUM, 0, xolotlComm);

			// If it is procId == 0 just store the integrated value in the
			// myPoints vector
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
		surfacePlotXY3D->getDataProvider()->setDataPoints(myPoints);

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
PetscErrorCode
monitorSurfaceXZ3D(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void*)
{
	// Initial declarations
	PetscErrorCode ierr;
	const double ****solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;
	double x, z;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % 10 != 0)
		PetscFunctionReturn(0);

	// Gets the process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the solver handler and local coordinates
	auto& solverHandler = PetscSolver::getSolverHandler();
	solverHandler.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the physical grid in the x direction
	auto grid = solverHandler.getXGrid();

	// Setup step size variables
	double hz = solverHandler.getStepSizeZ();

	// Choice of the cluster to be plotted
	IdType iCluster = 0;

	// Create a DataPoint vector to store the data to give to the data provider
	// for the visualization
	auto myPoints =
		std::make_shared<std::vector<viz::dataprovider::DataPoint>>();
	// Create a point here so that it is not created and deleted in the loop
	viz::dataprovider::DataPoint thePoint;

	// Loop on the full grid, Z and X first because they are the axis of the
	// plot
	for (auto k = 0; k < Mz; k++) {
		// Compute z
		z = (double)k * hz;

		for (auto i = 0; i < Mx; i++) {
			// Compute x
			x = (grid[i] + grid[i + 1]) / 2.0 - grid[1];

			// Initialize the value of the concentration to integrate over Y
			double conc = 0.0;

			for (auto j = 0; j < My; j++) {
				// If it is the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym &&
					k >= zs && k < zs + zm) {
					// Get the pointer to the beginning of the solution data for
					// this grid point
					gridPointSolution = solutionArray[k][j][i];

					// Integrate over Y
					conc += gridPointSolution[iCluster];
				}
			} // End of the loop on Y

			// Sum all the concentration on Y
			double totalConc = 0.0;
			MPI_Reduce(
				&conc, &totalConc, 1, MPI_DOUBLE, MPI_SUM, 0, xolotlComm);

			// If it is procId == 0 just store the integrated value in the
			// myPoints vector
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
		surfacePlotXZ3D->getDataProvider()->setDataPoints(myPoints);

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
PetscErrorCode
eventFunction3D(TS ts, PetscReal time, Vec solution, PetscScalar* fvalue, void*)
{
	// Initial declaration
	PetscErrorCode ierr;
	double ****solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;
	fvalue[0] = 1.0, fvalue[1] = 1.0;
	depthPositions3D.clear();

	PetscFunctionBeginUser;

	PetscInt TSNumber = -1;
	ierr = TSGetStepNumber(ts, &TSNumber);

	// Skip if it is the same TS as before
	if (TSNumber == previousTSNumber3D)
		PetscFunctionReturn(0);

	// Set the previous TS number
	previousTSNumber3D = TSNumber;

	// Gets the process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the solver handler and local coordinates
	auto& solverHandler = PetscSolver::getSolverHandler();
	solverHandler.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the network
	using NetworkType = core::network::IPSIReactionNetwork;
	auto& network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	// Get the number of species
	auto numSpecies = network.getSpeciesListSize();
	auto specIdI = network.getInterstitialSpeciesId();

	// Get the physical grid and step size
	auto grid = solverHandler.getXGrid();
	double hy = solverHandler.getStepSizeY();
	double hz = solverHandler.getStepSizeZ();

	// Get the flux handler to know the flux amplitude.
	auto fluxHandler = solverHandler.getFluxHandler();
	double heliumFluxAmplitude = fluxHandler->getFluxAmplitude();

	// Get the delta time from the previous timestep to this timestep
	double dt = time - solverHandler.getPreviousTime();

	// Work of the moving surface first
	if (solverHandler.moveSurface()) {
		// Write the initial surface positions
		if (procId == 0 && util::equal(time, 0.0)) {
			std::ofstream outputFile;
			outputFile.open("surface.txt", std::ios::app);
			outputFile << time << " ";

			// Loop on the possible yj
			for (auto yj = 0; yj < My; yj++) {
				for (auto zk = 0; zk < Mz; zk++) {
					// Get the position of the surface at yj, zk
					auto surfacePos = solverHandler.getSurfacePosition(yj, zk);
					outputFile << (double)yj * hy << " " << (double)zk * hz
							   << " " << grid[surfacePos + 1] - grid[1] << " ";
				}
			}
			outputFile << std::endl;
			outputFile.close();
		}

		// Get the initial vacancy concentration
		double initialVConc = solverHandler.getInitialVConc();

		// Loop on the possible zk and yj
		for (auto zk = 0; zk < Mz; zk++) {
			for (auto yj = 0; yj < My; yj++) {
				// Compute the total density of intersitials that escaped from
				// the surface since last timestep using the stored flux
				nInterstitial3D[yj][zk] += previousIFlux3D[yj][zk] * dt;

				// Remove the sputtering yield since last timestep
				nInterstitial3D[yj][zk] -=
					sputteringYield3D * heliumFluxAmplitude * dt;

				// Get the position of the surface at yj
				auto surfacePos = solverHandler.getSurfacePosition(yj, zk);
				auto xi = surfacePos + solverHandler.getLeftOffset();

				// Initialize the value for the flux
				auto myFlux = std::vector<double>(numSpecies, 0.0);

				// if xi is on this process
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym &&
					zk >= zs && zk < zs + zm) {
					// Get the concentrations at xi = surfacePos + 1
					gridPointSolution = solutionArray[zk][yj][xi];

					// Factor for finite difference
					double hxLeft = 0.0, hxRight = 0.0;
					if (xi >= 1 && xi < Mx) {
						hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
						hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
					}
					else if (xi < 1) {
						hxLeft = grid[xi + 1] - grid[xi];
						hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
					}
					else {
						hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
						hxRight = grid[xi + 1] - grid[xi];
					}
					double factor = hy * hz * 2.0 / (hxLeft + hxRight);

					network.updateOutgoingDiffFluxes(gridPointSolution, factor,
						iClusterIds3D, myFlux, xi - xs);
				}

				// Gather newFlux values at this position
				double newTotalFlux = 0.0;
				MPI_Allreduce(&myFlux[specIdI()], &newTotalFlux, 1, MPI_DOUBLE,
					MPI_SUM, xolotlComm);

				// Update the previous flux
				previousIFlux3D[yj][zk] = newTotalFlux;

				// Compare nInterstitials to the threshold to know if we should
				// move the surface

				// The density of tungsten is 62.8 atoms/nm3, thus the threshold
				// is
				double threshold =
					(62.8 - initialVConc) * (grid[xi] - grid[xi - 1]) * hy * hz;
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
		using NetworkType = core::network::IPSIReactionNetwork;
		auto psiNetwork = dynamic_cast<NetworkType*>(&network);
		auto dof = network.getDOF();
		auto specIdHe = psiNetwork->getHeliumSpeciesId();

		// Compute the prefactor for the probability (arbitrary)
		double prefactor =
			heliumFluxAmplitude * dt * solverHandler.getBurstingFactor();

		// The depth parameter to know where the bursting should happen
		double depthParam = solverHandler.getTauBursting(); // nm
		// The number of He per V in a bubble
		double heVRatio = solverHandler.getHeVRatio();

		// For now we are not bursting
		bool burst = false;

		// Loop on the full grid
		for (auto zk = 0; zk < Mz; zk++) {
			for (auto yj = 0; yj < My; yj++) {
				// Get the surface position
				auto surfacePos = solverHandler.getSurfacePosition(yj, zk);
				for (auto xi = surfacePos + solverHandler.getLeftOffset();
					 xi < Mx - solverHandler.getRightOffset(); xi++) {
					// If this is the locally owned part of the grid
					if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym &&
						zk >= zs && zk < zs + zm) {
						// Get the pointer to the beginning of the solution data
						// for this grid point
						gridPointSolution = solutionArray[zk][yj][xi];

						using HostUnmanaged = Kokkos::View<double*,
							Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
						auto hConcs = HostUnmanaged(gridPointSolution, dof);
						auto dConcs =
							Kokkos::View<double*>("Concentrations", dof);
						deep_copy(dConcs, hConcs);

						// Get the distance from the surface
						double distance = (grid[xi] + grid[xi + 1]) / 2.0 -
							grid[surfacePos + 1];

						// Compute the helium density at this grid point
						double heDensity =
							psiNetwork->getTotalAtomConcentration(
								dConcs, specIdHe, 1);

						// Compute the radius of the bubble from the number of
						// helium
						double nV = heDensity * (grid[xi + 1] - grid[xi]) * hy *
							hz / heVRatio;
						double latticeParam = network.getLatticeParameter();
						double tlcCubed =
							latticeParam * latticeParam * latticeParam;
						double radius = (sqrt(3.0) / 4) * latticeParam +
							cbrt((3.0 * tlcCubed * nV) / (8.0 * core::pi)) -
							cbrt((3.0 * tlcCubed) / (8.0 * core::pi));

						// If the radius is larger than the distance to the
						// surface, burst
						if (radius > distance) {
							burst = true;
							depthPositions3D.push_back(
								{static_cast<PetscInt>(zk),
									static_cast<PetscInt>(yj),
									static_cast<PetscInt>(xi)});
							// Exit the loop
							continue;
						}
						// Add randomness
						double prob = prefactor *
							(1.0 - (distance - radius) / distance) *
							std::min(1.0,
								exp(-(distance - depthParam) /
									(depthParam * 2.0)));
						double test = solverHandler.getRNG().GetRandomDouble();

						if (prob > test) {
							burst = true;
							depthPositions3D.push_back(
								{static_cast<PetscInt>(zk),
									static_cast<PetscInt>(yj),
									static_cast<PetscInt>(xi)});
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
PetscErrorCode
postEventFunction3D(TS ts, PetscInt nevents, PetscInt eventList[],
	PetscReal time, Vec solution, PetscBool, void*)
{
	// Initial declaration
	PetscErrorCode ierr;
	double ****solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Check if the surface has moved
	if (nevents == 0)
		PetscFunctionReturn(0);

	// Gets the process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the solver handler and local coordinates
	auto& solverHandler = PetscSolver::getSolverHandler();
	solverHandler.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the network
	auto& network = solverHandler.getNetwork();
	auto dof = network.getDOF();

	// Get the physical grid
	auto grid = solverHandler.getXGrid();
	// Get the step sizes
	double hy = solverHandler.getStepSizeY();
	double hz = solverHandler.getStepSizeZ();

	// Take care of bursting
	using NetworkType = core::network::IPSIReactionNetwork;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);

	// Loop on each bursting depth
	for (auto i = 0; i < depthPositions3D.size(); i++) {
		// Get the coordinates of the point
		auto xi = std::get<2>(depthPositions3D[i]),
			 yj = std::get<1>(depthPositions3D[i]),
			 zk = std::get<0>(depthPositions3D[i]);
		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[zk][yj][xi];

		// Get the surface position
		auto surfacePos = solverHandler.getSurfacePosition(yj, zk);
		// Get the distance from the surface
		double distance =
			(grid[xi] + grid[xi + 1]) / 2.0 - grid[surfacePos + 1];

		std::cout << "bursting at: " << zk * hz << " " << yj * hy << " "
				  << distance << std::endl;

		// Pinhole case
		auto nBurst = std::vector<double>(3, 0.0); // Not actually used here
		psiNetwork->updateBurstingConcs(gridPointSolution, 0.0, nBurst);
	}

	// Now takes care of moving surface
	bool moving = false;
	for (auto i = 0; i < nevents; i++) {
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
	for (auto zk = 0; zk < Mz; zk++) {
		for (auto yj = 0; yj < My; yj++) {
			// Get the position of the surface at yj
			auto surfacePos = solverHandler.getSurfacePosition(yj, zk);
			auto xi = surfacePos + solverHandler.getLeftOffset();

			// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
			double threshold =
				(62.8 - initialVConc) * (grid[xi] - grid[xi - 1]) * hy * hz;

			// Move the surface up
			if (nInterstitial3D[yj][zk] > threshold) {
				int nGridPoints = 0;
				// Move the surface up until it is smaller than the next
				// threshold
				while (nInterstitial3D[yj][zk] > threshold &&
					surfacePos + solverHandler.getLeftOffset() - 2 >= 0) {
					// Move the surface higher
					surfacePos--;
					xi = surfacePos + solverHandler.getLeftOffset();
					nGridPoints++;
					// Update the number of interstitials
					nInterstitial3D[yj][zk] -= threshold;
					// Update the thresold
					threshold = (62.8 - initialVConc) *
						(grid[xi] - grid[xi - 1]) * hy * hz;
				}

				// Throw an exception if the position is negative
				if (surfacePos + solverHandler.getLeftOffset() < 2) {
					PetscBool flagCheck;
					ierr = PetscOptionsHasName(
						NULL, NULL, "-check_collapse", &flagCheck);
					CHKERRQ(ierr);
					if (flagCheck) {
						// Write the convergence reason
						std::ofstream outputFile;
						outputFile.open("solverStatus.txt");
						outputFile << "overgrid" << std::endl;
						outputFile.close();
					}
					throw std::runtime_error(
						"\nxolotlSolver::Monitor3D: The surface "
						"is trying to go outside of the grid!!");
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

				// Initialize the vacancy concentration and the temperature on
				// the new grid points Get the single vacancy ID
				auto singleVacancyCluster = network.getSingleVacancy();
				auto vacancyIndex =
					core::network::IReactionNetwork::invalidIndex();
				if (singleVacancyCluster.getId() !=
					core::network::IReactionNetwork::invalidIndex())
					vacancyIndex = singleVacancyCluster.getId();
				// Get the surface temperature
				double temp = 0.0;
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym &&
					zk >= zs && zk < zs + zm) {
					temp = solutionArray[zk][yj][xi][dof];
				}
				double surfTemp = 0.0;
				MPI_Allreduce(
					&temp, &surfTemp, 1, MPI_DOUBLE, MPI_SUM, xolotlComm);
				// Loop on the new grid points
				while (nGridPoints >= 0) {
					// Position of the newly created grid point
					xi = surfacePos + nGridPoints;

					// If xi is on this process
					if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym &&
						zk >= zs && zk < zs + zm) {
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
				// Move it back as long as the number of interstitials in
				// negative
				while (nInterstitial3D[yj][zk] < 0.0) {
					// Compute the threshold to a deeper grid point
					threshold = (62.8 - initialVConc) *
						(grid[xi + 1] - grid[xi]) * hy * hz;
					// Set all the concentrations to 0.0 at xi = surfacePos + 1
					// if xi is on this process
					if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym &&
						zk >= zs && zk < zs + zm) {
						// Get the concentrations at xi = surfacePos + 1
						gridPointSolution = solutionArray[zk][yj][xi];
						// Loop on DOF
						for (auto i = 0; i < dof; i++) {
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
							  << yj * hy << " and " << zk * hz
							  << " at time: " << time << " s." << std::endl;
				}

				// Set it in the solver
				solverHandler.setSurfacePosition(surfacePos, yj, zk);
			}
		}
	}

	// Write the surface positions
	if (procId == 0) {
		std::ofstream outputFile;
		outputFile.open("surface.txt", std::ios::app);
		outputFile << time << " ";

		// Loop on the possible yj
		for (auto yj = 0; yj < My; yj++) {
			for (auto zk = 0; zk < Mz; zk++) {
				// Get the position of the surface at yj, zk
				auto surfacePos = solverHandler.getSurfacePosition(yj, zk);
				outputFile << (double)yj * hy << " " << (double)zk * hz << " "
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
PetscErrorCode
setupPetsc3DMonitor(TS ts)
{
	PetscErrorCode ierr;

	// Get the process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the xolotlViz handler registry
	auto vizHandlerRegistry = viz::VizHandlerRegistry::get();

	// Flags to launch the monitors or not
	PetscBool flagCheck, flagPerf, flagHeRetention, flagXeRetention, flagStatus,
		flag2DXYPlot, flag2DXZPlot, flagLargest;

	// Check the option -check_collapse
	ierr = PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCheck);
	checkPetscError(ierr,
		"setupPetsc3DMonitor: PetscOptionsHasName (-check_collapse) failed.");

	// Check the option -plot_perf
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_perf", &flagPerf);
	checkPetscError(
		ierr, "setupPetsc3DMonitor: PetscOptionsHasName (-plot_perf) failed.");

	// Check the option -plot_2d_xy
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_2d_xy", &flag2DXYPlot);
	checkPetscError(
		ierr, "setupPetsc3DMonitor: PetscOptionsHasName (-plot_2d_xy) failed.");

	// Check the option -plot_2d_xz
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_2d_xz", &flag2DXZPlot);
	checkPetscError(
		ierr, "setupPetsc3DMonitor: PetscOptionsHasName (-plot_2d_xz) failed.");

	// Check the option -helium_retention
	ierr =
		PetscOptionsHasName(NULL, NULL, "-helium_retention", &flagHeRetention);
	checkPetscError(ierr,
		"setupPetsc3DMonitor: PetscOptionsHasName (-helium_retention) failed.");

	// Check the option -xenon_retention
	ierr =
		PetscOptionsHasName(NULL, NULL, "-xenon_retention", &flagXeRetention);
	checkPetscError(ierr,
		"setupPetsc3DMonitor: PetscOptionsHasName (-xenon_retention) failed.");

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, NULL, "-start_stop", &flagStatus);
	checkPetscError(
		ierr, "setupPetsc3DMonitor: PetscOptionsHasName (-start_stop) failed.");

	// Check the option -largest_conc
	ierr = PetscOptionsHasName(NULL, NULL, "-largest_conc", &flagLargest);
	checkPetscError(ierr,
		"setupPetsc3DMonitor: PetscOptionsHasName (-largest_conc) failed.");

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto& network = solverHandler.getNetwork();

	// Determine if we have an existing restart file,
	// and if so, it it has had timesteps written to it.
	std::unique_ptr<io::XFile> networkFile;
	std::unique_ptr<io::XFile::TimestepGroup> lastTsGroup;

	std::string networkName = solverHandler.getNetworkName();
	bool hasConcentrations = false;
	if (not networkName.empty()) {
		networkFile = std::make_unique<io::XFile>(networkName);
		auto concGroup = networkFile->getGroup<io::XFile::ConcentrationGroup>();
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
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);
	checkPetscError(ierr, "setupPetsc3DMonitor: DMDAGetInfo failed.");

	// Set the post step processing to stop the solver if the time step
	// collapses
	if (flagCheck) {
		// Find the threshold
		PetscBool flag;
		ierr = PetscOptionsGetReal(
			NULL, NULL, "-check_collapse", &timeStepThreshold, &flag);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: PetscOptionsGetReal (-check_collapse) "
			"failed.");
		if (!flag)
			timeStepThreshold = 1.0e-16;

		// Set the post step process that tells the solver when to stop if the
		// time step collapse
		ierr = TSSetPostStep(ts, checkTimeStep);
		checkPetscError(
			ierr, "setupPetsc3DMonitor: TSSetPostStep (checkTimeStep) failed.");
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr = PetscOptionsGetReal(
			NULL, NULL, "-start_stop", &hdf5Stride3D, &flag);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: PetscOptionsGetReal (-start_stop) failed.");
		if (!flag)
			hdf5Stride3D = 1.0;

		// Compute the correct hdf5Previous3D for a restart
		if (hasConcentrations) {
			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			solverHandler.setPreviousTime(previousTime);
			hdf5Previous3D = (PetscInt)(previousTime / hdf5Stride3D);
		}

		// Don't do anything if both files have the same name
		if (hdf5OutputName3D != solverHandler.getNetworkName()) {
			// Get the solver handler
			auto& solverHandler = PetscSolver::getSolverHandler();

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
				io::XFile checkpointFile(
					hdf5OutputName3D, grid, xolotlComm, My, hy, Mz, hz);
			}

			// Copy the network group from the given file (if it has one).
			// We open the files using a single-process MPI communicator
			// because it is faster for a single process to do the
			// copy with HDF5's H5Ocopy implementation than it is
			// when all processes call the copy function.
			// The checkpoint file must be closed before doing this.
			writeNetwork(xolotlComm, solverHandler.getNetworkName(),
				hdf5OutputName3D, network);
		}

		// startStop3D will be called at each timestep
		ierr = TSMonitorSet(ts, startStop3D, NULL, NULL);
		checkPetscError(
			ierr, "setupPetsc3DMonitor: TSMonitorSet (startStop3D) failed.");
	}

	// If the user wants the surface to be able to move
	if (solverHandler.moveSurface() || solverHandler.burstBubbles()) {
		// Surface
		if (solverHandler.moveSurface()) {
			using NetworkType = core::network::IPSIReactionNetwork;
			using AmountType = NetworkType::AmountType;
			auto psiNetwork = dynamic_cast<NetworkType*>(&network);
			// Get the number of species
			auto numSpecies = psiNetwork->getSpeciesListSize();
			auto specIdI = psiNetwork->getInterstitialSpeciesId();

			// Initialize the composition
			auto comp = std::vector<AmountType>(numSpecies, 0);

			// Loop on interstital clusters
			bool iClusterExists = true;
			AmountType iSize = 1;
			while (iClusterExists) {
				comp[specIdI()] = iSize;
				auto clusterId = psiNetwork->findClusterId(comp);
				// Check that the helium cluster is present in the network
				if (clusterId != NetworkType::invalidIndex()) {
					iClusterIds3D.push_back(clusterId);
					iSize++;
				}
				else
					iClusterExists = false;
			}

			// Initialize nInterstitial3D and previousIFlux3D before monitoring
			// the interstitial flux
			for (auto j = 0; j < My; j++) {
				// Create a one dimensional vector of double
				std::vector<double> tempVector;
				for (auto k = 0; k < Mz; k++) {
					tempVector.push_back(0.0);
				}
				// Add the tempVector to nInterstitial3D and previousIFlux3D
				// to create their initial structure
				nInterstitial3D.push_back(tempVector);
				previousIFlux3D.push_back(tempVector);
			}

			// Get the interstitial information at the surface if concentrations
			// were stored
			if (hasConcentrations) {
				// Get the interstitial quantity from the HDF5 file
				nInterstitial3D = lastTsGroup->readData3D("nInterstitial");
				// Get the previous I flux from the HDF5 file
				previousIFlux3D = lastTsGroup->readData3D("previousFluxI");
				// Get the previous time from the HDF5 file
				double previousTime = lastTsGroup->readPreviousTime();
				solverHandler.setPreviousTime(previousTime);
			}

			// Get the sputtering yield
			sputteringYield3D = solverHandler.getSputteringYield();

			// Master process
			if (procId == 0) {
				// Clear the file where the surface will be written
				std::ofstream outputFile;
				outputFile.open("surface.txt");
				outputFile << "#time heights" << std::endl;
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
			perfPlot = vizHandlerRegistry->getPlot(viz::PlotType::SCATTER);

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>();
			labelProvider->axis1Label = "Process ID";
			labelProvider->axis2Label = "Solver Time";

			// Give it to the plot
			perfPlot->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider =
				std::make_shared<viz::dataprovider::CvsXDataProvider>();

			// Give it to the plot
			perfPlot->setDataProvider(dataProvider);
		}

		// monitorPerf will be called at each timestep
		ierr = TSMonitorSet(ts, monitorPerf, NULL, NULL);
		checkPetscError(
			ierr, "setupPetsc3DMonitor: TSMonitorSet (monitorPerf) failed.");
	}

	// Set the monitor to compute the helium fluence for the retention
	// calculation
	if (flagHeRetention) {
		// Get the previous time if concentrations were stored and initialize
		// the fluence
		if (hasConcentrations) {
			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			solverHandler.setPreviousTime(previousTime);
			// Initialize the fluence
			auto fluxHandler = solverHandler.getFluxHandler();
			// Increment the fluence with the value at this current timestep
			fluxHandler->computeFluence(previousTime);
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeFluence, NULL, NULL);
		checkPetscError(
			ierr, "setupPetsc3DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeHeliumRetention3D will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumRetention3D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: TSMonitorSet (computeHeliumRetention3D) "
			"failed.");

		// Master process
		if (procId == 0) {
			auto numSpecies = network.getSpeciesListSize();
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("retentionOut.txt");
			outputFile << "#fluence ";
			for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
				auto speciesName = network.getSpeciesName(id);
				outputFile << speciesName << "_content ";
			}
			outputFile << std::endl;
			outputFile.close();
		}
	}

	// Set the monitor to compute the xenon fluence and the retention
	// for the retention calculation
	if (flagXeRetention) {
		// Get the da from ts
		DM da;
		ierr = TSGetDM(ts, &da);
		checkPetscError(ierr, "setupPetsc3DMonitor: TSGetDM failed.");
		// Get the local boundaries
		PetscInt xm, ym, zm;
		ierr = DMDAGetCorners(da, NULL, NULL, NULL, &xm, &ym, &zm);
		checkPetscError(ierr, "setupPetsc3DMonitor: DMDAGetCorners failed.");
		// Create the local vectors on each process
		solverHandler.createLocalNE(xm, ym, zm);

		// Get the previous time if concentrations were stored and initialize
		// the fluence
		if (hasConcentrations) {
			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			solverHandler.setPreviousTime(previousTime);
			// Initialize the fluence
			auto fluxHandler = solverHandler.getFluxHandler();
			// Increment the fluence with the value at this current timestep
			fluxHandler->computeFluence(previousTime);
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeFluence, NULL, NULL);
		checkPetscError(
			ierr, "setupPetsc3DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeXenonRetention3D will be called at each timestep
		ierr = TSMonitorSet(ts, computeXenonRetention3D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: TSMonitorSet (computeXenonRetention3D) "
			"failed.");

		// Master process
		if (procId == 0) {
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("retentionOut.txt");
			outputFile << "#time Xenon_content radius partial_radius "
						  "Xenon_gb"
					   << std::endl;
			outputFile.close();
		}
	}

	// Set the monitor to save surface plots of clusters concentration
	if (flag2DXYPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a SurfacePlot
			surfacePlotXY3D =
				vizHandlerRegistry->getPlot(viz::PlotType::SURFACE);

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>();
			labelProvider->axis1Label = "Depth (nm)";
			labelProvider->axis2Label = "Y (nm)";
			labelProvider->axis3Label = "Concentration";

			// Give it to the plot
			surfacePlotXY3D->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider =
				std::make_shared<viz::dataprovider::CvsXYDataProvider>();

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
			surfacePlotXZ3D =
				vizHandlerRegistry->getPlot(viz::PlotType::SURFACE);

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>();
			labelProvider->axis1Label = "Depth (nm)";
			labelProvider->axis2Label = "Z (nm)";
			labelProvider->axis3Label = "Concentration";

			// Give it to the plot
			surfacePlotXZ3D->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider =
				std::make_shared<viz::dataprovider::CvsXYDataProvider>();

			// Give it to the plot
			surfacePlotXZ3D->setDataProvider(dataProvider);
		}

		// monitorSurfaceXZ3D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSurfaceXZ3D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: TSMonitorSet (monitorSurfaceXZ3D) failed.");
	}

	// Set the monitor to monitor the concentration of the largest cluster
	if (flagLargest) {
		// Look for the largest cluster
		auto& network = solverHandler.getNetwork();
		largestClusterId3D = network.getLargestClusterId();

		// Find the threshold
		PetscBool flag;
		ierr = PetscOptionsGetReal(
			NULL, NULL, "-largest_conc", &largestThreshold3D, &flag);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: PetscOptionsGetReal (-largest_conc) failed.");

		// monitorLargest3D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorLargest3D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: TSMonitorSet (monitorLargest3D) failed.");
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(ts, monitorTime, NULL, NULL);
	checkPetscError(
		ierr, "setupPetsc3DMonitor: TSMonitorSet (monitorTime) failed.");

	PetscFunctionReturn(0);
}

} /* end namespace monitor */
} /* end namespace solver */
} /* end namespace xolotl */
