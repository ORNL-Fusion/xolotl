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
#include <xolotl/core/network/AlloyReactionNetwork.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/io/XFile.h>
#include <xolotl/perf/xolotlPerf.h>
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

//! The pointer to the plot used in monitorScatter1D.
std::shared_ptr<viz::IPlot> scatterPlot1D;
//! The pointer to the series plot used in monitorSeries1D.
std::shared_ptr<viz::IPlot> seriesPlot1D;
//! The pointer to the 2D plot used in MonitorSurface.
std::shared_ptr<viz::IPlot> surfacePlot1D;
//! The variable to store the particle flux at the previous time step.
std::vector<double> previousSurfFlux1D, previousBulkFlux1D;
double previousIEventFlux1D = 0.0;
//! The variable to store the total number of atoms going through the surface or
//! bottom.
std::vector<double> nSurf1D, nBulk1D;
double nInterEvent1D = 0.0, nHeliumBurst1D = 0.0, nDeuteriumBurst1D = 0.0,
	   nTritiumBurst1D = 0.0;
//! The variable to store the xenon flux at the previous time step.
double previousXeFlux1D = 0.0;
//! The variable to store the sputtering yield at the surface.
double sputteringYield1D = 0.0;
//! The threshold for the negative concentration
double negThreshold1D = 0.0;
//! How often HDF5 file is written
PetscReal hdf5Stride1D = 0.0;
//! Previous time for HDF5
PetscInt hdf5Previous1D = 0;
//! HDF5 output file name
std::string hdf5OutputName1D = "xolotlStop.h5";
// The vector of depths at which bursting happens
std::vector<int> depthPositions1D;
// The vector of ids for diffusing interstitial clusters
std::vector<int> iClusterIds1D;
// Tracks the previous TS number
int previousTSNumber1D = -1;
// The id of the largest cluster
int largestClusterId1D = -1;
// The concentration threshold for the largest cluster
double largestThreshold1D = 1.0e-12;

// Timers
std::shared_ptr<perf::ITimer> initTimer;
std::shared_ptr<perf::ITimer> checkNegativeTimer;
std::shared_ptr<perf::ITimer> tridynTimer;
std::shared_ptr<perf::ITimer> startStopTimer;
std::shared_ptr<perf::ITimer> heRetentionTimer;
std::shared_ptr<perf::ITimer> xeRetentionTimer;
std::shared_ptr<perf::ITimer> scatterTimer;
std::shared_ptr<perf::ITimer> seriesTimer;
std::shared_ptr<perf::ITimer> eventFuncTimer;
std::shared_ptr<perf::ITimer> postEventFuncTimer;

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "checkNegative1D")
/**
 * This is a monitoring method that looks at if there are negative
 * concentrations at each time step.
 */
PetscErrorCode
checkNegative1D(TS ts, PetscInt timestep, PetscReal time, Vec solution, void*)
{
	perf::ScopedTimer myTimer(checkNegativeTimer);

	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the MPI communicator
	auto xolotlComm = util::getMPIComm();
	// Get the number of processes
	int worldSize;
	MPI_Comm_size(xolotlComm, &worldSize);
	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network and dof
	auto& network = solverHandler.getNetwork();
	const int nClusters = network.getNumClusters();

	// Loop on the local grid
	for (PetscInt i = xs; i < xs + xm; i++) {
		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[i]; // Loop on the concentrations
		for (int l = 0; l < nClusters; l++) {
			if (gridPointSolution[l] < negThreshold1D &&
				gridPointSolution[l] > 0.0) {
				gridPointSolution[l] = negThreshold1D;
			}
			else if (gridPointSolution[l] > -negThreshold1D &&
				gridPointSolution[l] < 0.0) {
				gridPointSolution[l] = -negThreshold1D;
			}
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorLargest1D")
/**
 * This is a monitoring method that looks at the largest cluster concentration
 */
PetscErrorCode
monitorLargest1D(TS ts, PetscInt timestep, PetscReal time, Vec solution, void*)
{
	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the MPI communicator
	auto xolotlComm = util::getMPIComm();
	// Get the number of processes
	int worldSize;
	MPI_Comm_size(xolotlComm, &worldSize);
	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Loop on the local grid
	for (PetscInt i = xs; i < xs + xm; i++) {
		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[i];
		// Check the concentration
		if (gridPointSolution[largestClusterId1D] > largestThreshold1D) {
			ierr = TSSetConvergedReason(ts, TS_CONVERGED_USER);
			CHKERRQ(ierr);
			// Send an error
			throw std::string("\nxolotlSolver::Monitor1D: The largest cluster "
							  "concentration is too high!!");
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeTRIDYN1D")
/**
 * This is a monitoring method that will compute the data to send to TRIDYN
 */
PetscErrorCode
computeTRIDYN1D(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	perf::ScopedTimer myTimer(tridynTimer);

	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the MPI communicator
	auto xolotlComm = util::getMPIComm();

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(xolotlComm, &worldSize);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	using NetworkType =
		core::network::PSIReactionNetwork<core::network::PSIFullSpeciesList>;
	using Spec = typename NetworkType::Species;
	auto& network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	int dof = network.getDOF();

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the total size of the grid
	PetscInt Mx;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// Get the array of concentration
	PetscReal** solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Save current concentrations as an HDF5 file.
	//
	// First create the file for parallel file access.
	std::ostringstream tdFileStr;
	tdFileStr << "TRIDYN_" << timestep << ".h5";
	io::HDF5File tdFile(tdFileStr.str(),
		io::HDF5File::AccessMode::CreateOrTruncateIfExists, xolotlComm, true);

	// Define a dataset for concentrations.
	// Everyone must create the dataset with the same shape.
	constexpr auto numConcSpecies = 5;
	constexpr auto numValsPerGridpoint = numConcSpecies + 2;
	const auto firstIdxToWrite = (surfacePos + solverHandler.getLeftOffset());
	const auto numGridpointsWithConcs = (Mx - firstIdxToWrite);
	io::HDF5File::SimpleDataSpace<2>::Dimensions concsDsetDims = {
		(hsize_t)numGridpointsWithConcs, numValsPerGridpoint};
	io::HDF5File::SimpleDataSpace<2> concsDsetSpace(concsDsetDims);

	const std::string concsDsetName = "concs";
	io::HDF5File::DataSet<double> concsDset(
		tdFile, concsDsetName, concsDsetSpace);

	// Specify the concentrations we will write.
	// We only consider our own grid points.
	const auto myFirstIdxToWrite = std::max(xs, firstIdxToWrite);
	auto myEndIdx = (xs + xm); // "end" in the C++ sense; i.e., one-past-last
	auto myNumPointsToWrite =
		(myEndIdx > myFirstIdxToWrite) ? (myEndIdx - myFirstIdxToWrite) : 0;
	io::HDF5File::DataSet<double>::DataType2D<numValsPerGridpoint> myConcs(
		myNumPointsToWrite);

	for (auto xi = myFirstIdxToWrite; xi < myEndIdx; ++xi) {
		if (xi >= firstIdxToWrite) {
			// Determine current gridpoint value.
			double x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];

			// Access the solution data for this grid point.
			auto gridPointSolution = solutionArray[xi];
			using HostUnmanaged = Kokkos::View<double*, Kokkos::HostSpace,
				Kokkos::MemoryUnmanaged>;
			auto hConcs = HostUnmanaged(gridPointSolution, dof);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof);
			deep_copy(dConcs, hConcs);

			// Get the total concentrations at this grid point
			auto currIdx = xi - myFirstIdxToWrite;
			myConcs[currIdx][0] = (x - (grid[surfacePos + 1] - grid[1]));
			myConcs[currIdx][1] =
				network.getTotalAtomConcentration(dConcs, Spec::He, 1);
			myConcs[currIdx][2] =
				network.getTotalAtomConcentration(dConcs, Spec::D, 1);
			myConcs[currIdx][3] =
				network.getTotalAtomConcentration(dConcs, Spec::T, 1);
			myConcs[currIdx][4] =
				network.getTotalAtomConcentration(dConcs, Spec::V, 1);
			myConcs[currIdx][5] =
				network.getTotalAtomConcentration(dConcs, Spec::I, 1);
			myConcs[currIdx][6] = gridPointSolution[dof];
		}
	}

	// Write the concs dataset in parallel.
	// (We write only our part.)
	concsDset.parWrite2D<numValsPerGridpoint>(
		xolotlComm, myFirstIdxToWrite - firstIdxToWrite, myConcs);

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "startStop1D")
/**
 * This is a monitoring method that update an hdf5 file at each time step.
 */
PetscErrorCode
startStop1D(TS ts, PetscInt timestep, PetscReal time, Vec solution, void*)
{
	perf::ScopedTimer myTimer(startStopTimer);

	// Initial declaration
	PetscErrorCode ierr;
	const double **solutionArray, *gridPointSolution;
	PetscInt xs, xm, Mx;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Compute the dt
	double previousTime = solverHandler.getPreviousTime();
	double dt = time - previousTime;

	// Don't do anything if it is not on the stride
	if (((int)((time + dt / 10.0) / hdf5Stride1D) <= hdf5Previous1D) &&
		timestep > 0) {
		PetscFunctionReturn(0);
	}

	// Update the previous time
	if ((int)((time + dt / 10.0) / hdf5Stride1D) > hdf5Previous1D)
		hdf5Previous1D++;

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

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the network and dof
	auto& network = solverHandler.getNetwork();
	const int dof = network.getDOF();

	// Create an array for the concentration
	double concArray[dof + 1][2];

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Open the existing HDF5 file
	io::XFile checkpointFile(
		hdf5OutputName1D, xolotlComm, io::XFile::AccessMode::OpenReadWrite);

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	CHKERRQ(ierr);

	// Add a concentration time step group for the current time step.
	auto concGroup = checkpointFile.getGroup<io::XFile::ConcentrationGroup>();
	assert(concGroup);
	auto tsGroup = concGroup->addTimestepGroup(
		timestep, time, previousTime, currentTimeStep);

	// Get the names of the species in the network
	auto numSpecies = network.getSpeciesListSize();
	std::vector<std::string> names;
	for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
		names.push_back(network.getSpeciesName(id));
	}

	if (solverHandler.moveSurface() || solverHandler.getLeftOffset() == 1) {
		// Write the surface positions and the associated interstitial
		// quantities in the concentration sub group
		tsGroup->writeSurface1D(surfacePos, nInterEvent1D, previousIEventFlux1D,
			nSurf1D, previousSurfFlux1D, names);
	}

	// Write the bottom impurity information if the bottom is a free surface
	if (solverHandler.getRightOffset() == 1)
		tsGroup->writeBottom1D(nBulk1D, previousBulkFlux1D, names);

	// Write the bursting information if the bubble bursting is used
	if (solverHandler.burstBubbles())
		tsGroup->writeBursting1D(
			nHeliumBurst1D, nDeuteriumBurst1D, nTritiumBurst1D);

	// Determine the concentration values we will write.
	// We only examine and collect the grid points we own.
	// TODO measure impact of us building the flattened representation
	// rather than a ragged 2D representation.
	io::XFile::TimestepGroup::Concs1DType concs(xm);
	for (auto i = 0; i < xm; ++i) {
		// Access the solution data for the current grid point.
		auto gridPointSolution = solutionArray[xs + i];

		for (auto l = 0; l < dof + 1; ++l) {
			if (std::fabs(gridPointSolution[l]) > 1.0e-16) {
				concs[i].emplace_back(l, gridPointSolution[l]);
			}
		}
	}

	// Write our concentration data to the current timestep group
	// in the HDF5 file.
	// We only write the data for the grid points we own.
	tsGroup->writeConcentrations(checkpointFile, xs, concs);

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	ierr = computeTRIDYN1D(ts, timestep, time, solution, NULL);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeHeliumRetention1D")
/**
 * This is a monitoring method that will compute the helium retention
 */
PetscErrorCode
computeHeliumRetention1D(TS ts, PetscInt, PetscReal time, Vec solution, void*)
{
	perf::ScopedTimer myTimer(heRetentionTimer);

	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the flux handler that will be used to know the fluence
	auto fluxHandler = solverHandler.getFluxHandler();
	// Get the diffusion handler
	auto diffusionHandler = solverHandler.getDiffusionHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the total size of the grid
	PetscInt Mx;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();
	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Get the network
	using NetworkType = core::network::IPSIReactionNetwork;
	using AmountType = NetworkType::AmountType;
	auto& network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	const int dof = network.getDOF();

	// Get the array of concentration
	PetscReal** solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration over the grid
	auto numSpecies = network.getSpeciesListSize();
	auto specIdHe = network.getHeliumSpeciesId();
	auto myConcData = std::vector<double>(numSpecies, 0.0);

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal* gridPointSolution;

	// Loop on the grid
	for (PetscInt xi = xs; xi < xs + xm; xi++) {
		// Boundary conditions
		if (xi < surfacePos + solverHandler.getLeftOffset() ||
			xi >= Mx - solverHandler.getRightOffset())
			continue;

		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[xi];

		double hx = grid[xi + 1] - grid[xi];

		using HostUnmanaged =
			Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
		auto hConcs = HostUnmanaged(gridPointSolution, dof);
		auto dConcs = Kokkos::View<double*>("Concentrations", dof);
		deep_copy(dConcs, hConcs);

		// Get the total concentrations at this grid point
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			myConcData[id()] +=
				network.getTotalAtomConcentration(dConcs, id, 1) * hx;
		}
	}

	// Get the current process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Determine total concentrations for He, D, T.
	auto totalConcData = std::vector<double>(numSpecies, 0.0);

	MPI_Reduce(myConcData.data(), totalConcData.data(), numSpecies, MPI_DOUBLE,
		MPI_SUM, 0, xolotlComm);

	// Get the delta time from the previous timestep to this timestep
	double previousTime = solverHandler.getPreviousTime();
	double dt = time - previousTime;

	// Look at the fluxes leaving the free surface
	if (solverHandler.getLeftOffset() == 1) {
		// Set the surface position
		int xi = surfacePos + 1;

		// Value to know on which processor is the surface
		int surfaceProc = 0;

		// Check we are on the right proc
		if (xi >= xs && xi < xs + xm) {
			// Compute the total number of impurities that left at the surface
			for (std::size_t i = 0; i < numSpecies; ++i) {
				nSurf1D[i] += previousSurfFlux1D[i] * dt;
			}
			auto myFluxData = std::vector<double>(numSpecies, 0.0);

			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[xi];

			// Factor for finite difference
			double hxLeft = 0.0, hxRight = 0.0;
			if (xi - 1 >= 0 && xi < Mx) {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			}
			else if (xi - 1 < 0) {
				hxLeft = grid[xi + 1] - grid[xi];
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			}
			else {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = grid[xi + 1] - grid[xi];
			}
			double factor = 2.0 / (hxLeft + hxRight);

			// Get the vector of diffusing clusters
			auto diffusingIds = diffusionHandler->getDiffusingIds();

			network.updateOutgoingDiffFluxes(
				gridPointSolution, factor, diffusingIds, myFluxData, xi - xs);

			// Take into account the surface advection
			// Get the surface advection handler
			auto advecHandler = solverHandler.getAdvectionHandler();
			// Get the sink strengths and advecting clusters
			auto sinkStrengths = advecHandler->getSinkStrengths();
			auto advecClusters = advecHandler->getAdvectingClusters();
			// Set the distance from the surface
			double distance = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1] -
				advecHandler->getLocation();

			network.updateOutgoingAdvecFluxes(gridPointSolution,
				3.0 /
					(core::kBoltzmann * distance * distance * distance *
						distance),
				advecClusters, sinkStrengths, myFluxData, xi - xs);

			for (std::size_t i = 0; i < numSpecies; ++i) {
				previousSurfFlux1D[i] = myFluxData[i];
			}

			// Set the surface processor
			surfaceProc = procId;
		}

		// Get which processor will send the information
		// TODO do we need to do this allreduce just to figure out
		// who owns the data?
		// And is it supposed to be a sum?   Why not a min?
		int surfaceId = 0;
		MPI_Allreduce(
			&surfaceProc, &surfaceId, 1, MPI_INT, MPI_SUM, xolotlComm);

		// Send the information about impurities
		// to the other processes
		std::vector<double> countFluxData;
		for (std::size_t i = 0; i < numSpecies; ++i) {
			countFluxData.push_back(nSurf1D[i]);
			countFluxData.push_back(previousSurfFlux1D[i]);
		}
		MPI_Bcast(countFluxData.data(), countFluxData.size(), MPI_DOUBLE,
			surfaceId, xolotlComm);

		// Extract impurity data from broadcast buffer.
		for (std::size_t i = 0; i < numSpecies; ++i) {
			nSurf1D[i] = countFluxData[2 * i];
			previousSurfFlux1D[i] = countFluxData[(2 * i) + 1];
		}
	}

	// Look at the fluxes going in the bulk if the bottom is a free surface
	if (solverHandler.getRightOffset() == 1) {
		// Set the bottom surface position
		int xi = Mx - 2;

		// Value to know on which processor is the bottom
		int bottomProc = 0;

		// Check we are on the right proc
		if (xi >= xs && xi < xs + xm) {
			// Compute the total number of impurities that went in the bulk
			for (std::size_t i = 0; i < numSpecies; ++i) {
				nBulk1D[i] += previousBulkFlux1D[i] * dt;
			}
			auto myFluxData = std::vector<double>(numSpecies, 0.0);

			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[xi];

			// Factor for finite difference
			double hxLeft = 0.0, hxRight = 0.0;
			if (xi - 1 >= 0 && xi < Mx) {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			}
			else if (xi - 1 < 0) {
				hxLeft = grid[xi + 1] - grid[xi];
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			}
			else {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = grid[xi + 1] - grid[xi];
			}
			double factor = 2.0 / (hxLeft + hxRight);

			// Get the vector of diffusing clusters
			auto diffusingIds = diffusionHandler->getDiffusingIds();

			network.updateOutgoingDiffFluxes(
				gridPointSolution, factor, diffusingIds, myFluxData, xi - xs);

			for (std::size_t i = 0; i < numSpecies; ++i) {
				previousBulkFlux1D[i] = myFluxData[i];
			}

			// Set the bottom processor
			bottomProc = procId;
		}

		// Get which processor will send the information
		// TODO do we need to do this allreduce just to figure out
		// who owns the data?
		// And is it supposed to be a sum?   Why not a min?
		int bottomId = 0;
		MPI_Allreduce(&bottomProc, &bottomId, 1, MPI_INT, MPI_SUM, xolotlComm);

		// Send the information about impurities
		// to the other processes
		std::vector<double> countFluxData;
		for (std::size_t i = 0; i < numSpecies; ++i) {
			countFluxData.push_back(nBulk1D[i]);
			countFluxData.push_back(previousBulkFlux1D[i]);
		}
		MPI_Bcast(countFluxData.data(), countFluxData.size(), MPI_DOUBLE,
			bottomId, xolotlComm);

		// Extract impurity data from broadcast buffer.
		for (std::size_t i = 0; i < numSpecies; ++i) {
			nBulk1D[i] = countFluxData[2 * i];
			previousBulkFlux1D[i] = countFluxData[(2 * i) + 1];
		}
	}

	// Master process
	if (procId == 0) {
		// Get the fluence
		double fluence = fluxHandler->getFluence();

		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			std::cout << network.getSpeciesName(id)
					  << " content = " << totalConcData[id()] << '\n';
		}
		std::cout << "Fluence = " << fluence << '\n' << std::endl;

		// Uncomment to write the retention and the fluence in a file
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt", std::ios::app);
		outputFile << fluence << ' ';
		for (std::size_t i = 0; i < numSpecies; ++i) {
			outputFile << totalConcData[i] << ' ';
		}
		if (solverHandler.getRightOffset() == 1) {
			for (std::size_t i = 0; i < numSpecies; ++i) {
				outputFile << nBulk1D[i] << ' ';
			}
		}
		if (solverHandler.getLeftOffset() == 1) {
			for (std::size_t i = 0; i < numSpecies; ++i) {
				outputFile << nSurf1D[i] << ' ';
			}
		}
		outputFile << nHeliumBurst1D << " " << nDeuteriumBurst1D << " "
				   << nTritiumBurst1D << std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeXenonRetention1D")
/**
 * This is a monitoring method that will compute the xenon retention
 */
PetscErrorCode
computeXenonRetention1D(TS ts, PetscInt, PetscReal time, Vec solution, void*)
{
	perf::ScopedTimer myTimer(xeRetentionTimer);

	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the total size of the grid
	PetscInt Mx;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	using NetworkType = core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Degrees of freedom is the total number of clusters in the network
	auto& network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	const int dof = network.getDOF();

	// Get the complete data array, including ghost cells
	Vec localSolution;
	ierr = DMGetLocalVector(da, &localSolution);
	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da, solution, INSERT_VALUES, localSolution);
	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, solution, INSERT_VALUES, localSolution);
	CHKERRQ(ierr);
	// Get the array of concentration
	PetscReal** solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, localSolution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration and other values over the grid
	double xeConcentration = 0.0, bubbleConcentration = 0.0, radii = 0.0,
		   partialBubbleConcentration = 0.0, partialRadii = 0.0,
		   partialSize = 0.0;

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal* gridPointSolution;

	// Get the minimum size for the radius
	auto minSizes = solverHandler.getMinSizes();

	// Get Xe_1
	Composition xeComp = Composition::zero();
	xeComp[Spec::Xe] = 1;
	auto xeCluster = network.findCluster(xeComp, plsm::onHost);
	auto xeId = xeCluster.getId();

	// Loop on the grid
	for (PetscInt xi = xs; xi < xs + xm; xi++) {
		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[xi];

		using HostUnmanaged =
			Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
		auto hConcs = HostUnmanaged(gridPointSolution, dof);
		auto dConcs = Kokkos::View<double*>("Concentrations", dof);
		deep_copy(dConcs, hConcs);

		// Initialize the volume fraction and hx
		double hx = grid[xi + 1] - grid[xi];

		// Get the concentrations
		xeConcentration +=
			network.getTotalAtomConcentration(dConcs, Spec::Xe, 1) * hx;
		bubbleConcentration +=
			network.getTotalConcentration(dConcs, Spec::Xe, 1) * hx;
		radii += network.getTotalRadiusConcentration(dConcs, Spec::Xe, 1) * hx;
		partialBubbleConcentration =
			network.getTotalConcentration(dConcs, Spec::Xe, minSizes[0]) * hx;
		partialRadii +=
			network.getTotalRadiusConcentration(dConcs, Spec::Xe, minSizes[0]) *
			hx;
		partialSize +=
			network.getTotalAtomConcentration(dConcs, Spec::Xe, minSizes[0]) *
			hx;

		// Set the volume fraction
		double volumeFrac =
			network.getTotalVolumeFraction(dConcs, Spec::Xe, minSizes[0]);
		solverHandler.setVolumeFraction(volumeFrac, xi - xs);
		// Set the monomer concentration
		solverHandler.setMonomerConc(
			gridPointSolution[xeCluster.getId()], xi - xs);
	}

	// Get the current process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Sum all the concentrations through MPI reduce
	std::array<double, 6> myConcData{xeConcentration, bubbleConcentration,
		radii, partialBubbleConcentration, partialRadii, partialSize};
	std::array<double, 6> totalConcData{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
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
		int xi = std::get<0>(pair);
		// Check we are on the right proc
		if (xi >= xs && xi < xs + xm) {
			double previousXeFlux = std::get<1>(localNE[xi - xs][0][0]);
			globalXeFlux += previousXeFlux * (grid[xi + 1] - grid[xi]);
			// Set the amount in the vector we keep
			solverHandler.setLocalXeRate(previousXeFlux * dt, xi - xs);
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
		int xi = std::get<0>(pair);

		// Factor for finite difference
		double hxLeft = 0.0, hxRight = 0.0;
		if (xi - 1 >= 0 && xi < Mx) {
			hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
			hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
		}
		else if (xi - 1 < 0) {
			hxLeft = grid[xi + 1] - grid[xi];
			hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
		}
		else {
			hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
			hxRight = grid[xi + 1] - grid[xi];
		}
		double factor = 2.0 / (hxLeft + hxRight);
		// Check we are on the right proc
		if (xi >= xs && xi < xs + xm) {
			// Left
			xi = std::get<0>(pair) - 1;
			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[xi];
			// Compute the flux coming from the left
			network.updateOutgoingDiffFluxes(gridPointSolution, factor / hxLeft,
				diffusingIds, myRate, xi + 1 - xs);

			// Right
			xi = std::get<0>(pair) + 1;
			gridPointSolution = solutionArray[xi];
			// Compute the flux coming from the right
			network.updateOutgoingDiffFluxes(gridPointSolution,
				factor / hxRight, diffusingIds, myRate, xi + 1 - xs);

			// Middle
			xi = std::get<0>(pair);
			solverHandler.setPreviousXeFlux(myRate[0], xi - xs);
		}
	}

	// Master process
	if (procId == 0) {
		// Get the number of xenon that went to the GB
		double nXenon = solverHandler.getNXeGB();

		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		std::cout << "Xenon concentration = " << totalConcData[0] << std::endl;
		std::cout << "Xenon GB = " << nXenon << std::endl << std::endl;

		// Make sure the average partial radius makes sense
		double averagePartialRadius = 0.0, averagePartialSize = 0.0;
		if (totalConcData[3] > 1.e-16) {
			averagePartialRadius = totalConcData[4] / totalConcData[3];
			averagePartialSize = totalConcData[5] / totalConcData[3];
		}

		// Uncomment to write the content in a file
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt", std::ios::app);
		outputFile << time << " " << totalConcData[0] << " "
				   << totalConcData[2] / totalConcData[1] << " "
				   << averagePartialRadius << " " << totalConcData[3] << " "
				   << averagePartialSize << std::endl;
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
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "profileTemperature1D")
/**
 * This is a monitoring method that will store the temperature profile
 */
PetscErrorCode
profileTemperature1D(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Gets the process ID (important when it is running in parallel)
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network and dof
	auto& network = solverHandler.getNetwork();
	const int dof = network.getDOF();

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the total size of the grid
	PetscInt Mx;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// Get the array of concentration
	PetscReal** solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal* gridPointSolution;

	// Create the output file
	std::ofstream outputFile;
	if (procId == 0) {
		outputFile.open("tempProf.txt", std::ios::app);
		outputFile << time;
	}

	// Loop on the entire grid
	for (int xi = surfacePos + solverHandler.getLeftOffset();
		 xi < Mx - solverHandler.getRightOffset(); xi++) {
		// Set x
		double x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];

		double localTemp = 0.0;
		// Check if this process is in charge of xi
		if (xi >= xs && xi < xs + xm) {
			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[xi];

			// Get the local temperature
			localTemp = gridPointSolution[dof];
		}

		// Get the value on procId = 0
		double temperature = 0.0;
		MPI_Reduce(
			&localTemp, &temperature, 1, MPI_DOUBLE, MPI_SUM, 0, xolotlComm);

		// The master process writes in the file
		if (procId == 0) {
			outputFile << " " << temperature;
		}
	}

	// Close the file
	if (procId == 0) {
		outputFile << std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeAlloy1D")
/**
 * This is a monitoring method that will compute average density and diameter
 */
PetscErrorCode
computeAlloy1D(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	// Initial declarations
	PetscErrorCode ierr;
	PetscInt xs, xm;

	PetscFunctionBeginUser;

	// Get the number of processes
	auto xolotlComm = util::getMPIComm();
	int worldSize;
	MPI_Comm_size(xolotlComm, &worldSize);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the physical grid and its length
	auto grid = solverHandler.getXGrid();
	int xSize = grid.size();

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the total size of the grid
	PetscInt Mx;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the array of concentration
	PetscReal** solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	using NetworkType = core::network::AlloyReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Degrees of freedom is the total number of clusters in the network
	auto& network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	const int dof = network.getDOF();

	// Initial declarations for the density and diameter
	double iDensity = 0.0, vDensity = 0.0, voidDensity = 0.0,
		   frankDensity = 0.0, faultedDensity = 0.0, perfectDensity = 0.0,
		   voidPartialDensity = 0.0, frankPartialDensity = 0.0,
		   faultedPartialDensity = 0.0, perfectPartialDensity = 0.0,
		   iDiameter = 0.0, vDiameter = 0.0, voidDiameter = 0.0,
		   frankDiameter = 0.0, faultedDiameter = 0.0, perfectDiameter = 0.0,
		   voidPartialDiameter = 0.0, frankPartialDiameter = 0.0,
		   faultedPartialDiameter = 0.0, perfectPartialDiameter = 0.0;

	// Get the minimum size for the loop densities and diameters
	auto minSizes = solverHandler.getMinSizes();

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal* gridPointSolution;

	// Loop on the grid
	for (PetscInt xi = xs; xi < xs + xm; xi++) {
		// Boundary conditions
		if (xi < surfacePos + solverHandler.getLeftOffset() ||
			xi == Mx - solverHandler.getRightOffset())
			continue;

		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[xi];

		using HostUnmanaged =
			Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
		auto hConcs = HostUnmanaged(gridPointSolution, dof);
		auto dConcs = Kokkos::View<double*>("Concentrations", dof);
		deep_copy(dConcs, hConcs);

		// I
		iDensity += network.getTotalConcentration(dConcs, Spec::I, 1);
		iDiameter +=
			2.0 * network.getTotalRadiusConcentration(dConcs, Spec::I, 1);

		// V
		vDensity += network.getTotalConcentration(dConcs, Spec::V, 1);
		vDiameter +=
			2.0 * network.getTotalRadiusConcentration(dConcs, Spec::V, 1);

		// Void
		voidDensity += network.getTotalConcentration(dConcs, Spec::Void, 1);
		voidDiameter +=
			2.0 * network.getTotalRadiusConcentration(dConcs, Spec::Void, 1);
		voidPartialDensity +=
			network.getTotalConcentration(dConcs, Spec::Void, minSizes[0]);
		voidPartialDiameter += 2.0 *
			network.getTotalRadiusConcentration(
				dConcs, Spec::Void, minSizes[0]);

		// Faulted
		faultedDensity +=
			network.getTotalConcentration(dConcs, Spec::Faulted, 1);
		faultedDiameter +=
			2.0 * network.getTotalRadiusConcentration(dConcs, Spec::Faulted, 1);
		faultedPartialDensity +=
			network.getTotalConcentration(dConcs, Spec::Faulted, minSizes[1]);
		faultedPartialDiameter += 2.0 *
			network.getTotalRadiusConcentration(
				dConcs, Spec::Faulted, minSizes[1]);

		// Perfect
		perfectDensity +=
			network.getTotalConcentration(dConcs, Spec::Perfect, 1);
		perfectDiameter +=
			2.0 * network.getTotalRadiusConcentration(dConcs, Spec::Perfect, 1);
		perfectPartialDensity +=
			network.getTotalConcentration(dConcs, Spec::Perfect, minSizes[2]);
		perfectPartialDiameter += 2.0 *
			network.getTotalRadiusConcentration(
				dConcs, Spec::Perfect, minSizes[2]);

		// Frank
		frankDensity += network.getTotalConcentration(dConcs, Spec::Frank, 1);
		frankDiameter +=
			2.0 * network.getTotalRadiusConcentration(dConcs, Spec::Frank, 1);
		frankPartialDensity +=
			network.getTotalConcentration(dConcs, Spec::Frank, minSizes[3]);
		frankPartialDiameter += 2.0 *
			network.getTotalRadiusConcentration(
				dConcs, Spec::Frank, minSizes[3]);
	}

	// Sum all the concentrations through MPI reduce
	std::array<double, 20> myConcData{iDensity, vDensity, voidDensity,
		frankDensity, faultedDensity, perfectDensity, voidPartialDensity,
		frankPartialDensity, faultedPartialDensity, perfectPartialDensity,
		iDiameter, vDiameter, voidDiameter, frankDiameter, faultedDiameter,
		perfectDiameter, voidPartialDiameter, frankPartialDiameter,
		faultedPartialDiameter, perfectPartialDiameter};
	std::array<double, 20> totalConcData{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	MPI_Reduce(myConcData.data(), totalConcData.data(), myConcData.size(),
		MPI_DOUBLE, MPI_SUM, 0, xolotlComm);

	// Average the data
	if (procId == 0) {
		double iTotalDensity = totalConcData[0],
			   vTotalDensity = totalConcData[1],
			   voidTotalDensity = totalConcData[2],
			   frankTotalDensity = totalConcData[3],
			   faultedTotalDensity = totalConcData[4],
			   perfectTotalDensity = totalConcData[5],
			   voidPartialTotalDensity = totalConcData[6],
			   frankPartialTotalDensity = totalConcData[7],
			   faultedPartialTotalDensity = totalConcData[8],
			   perfectPartialTotalDensity = totalConcData[9],
			   iTotalDiameter = totalConcData[10],
			   vTotalDiameter = totalConcData[11],
			   voidTotalDiameter = totalConcData[12],
			   frankTotalDiameter = totalConcData[13],
			   faultedTotalDiameter = totalConcData[14],
			   perfectTotalDiameter = totalConcData[13],
			   voidPartialTotalDiameter = totalConcData[16],
			   frankPartialTotalDiameter = totalConcData[17],
			   faultedPartialTotalDiameter = totalConcData[18],
			   perfectPartialTotalDiameter = totalConcData[19];
		iTotalDensity = iTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
		vTotalDensity = vTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
		voidTotalDensity = voidTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
		perfectTotalDensity =
			perfectTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
		faultedTotalDensity =
			faultedTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
		frankTotalDensity =
			frankTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
		voidPartialTotalDensity =
			voidPartialTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
		perfectPartialTotalDensity =
			perfectPartialTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
		faultedPartialTotalDensity =
			faultedPartialTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
		frankPartialTotalDensity =
			frankPartialTotalDensity / (grid[Mx] - grid[surfacePos + 1]);
		iTotalDiameter = iTotalDiameter /
			(iTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
		vTotalDiameter = vTotalDiameter /
			(vTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
		voidTotalDiameter = voidTotalDiameter /
			(voidTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
		perfectTotalDiameter = perfectTotalDiameter /
			(perfectTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
		faultedTotalDiameter = faultedTotalDiameter /
			(faultedTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
		frankTotalDiameter = frankTotalDiameter /
			(frankTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
		voidPartialTotalDiameter = voidPartialTotalDiameter /
			(voidPartialTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
		perfectPartialTotalDiameter = perfectPartialTotalDiameter /
			(perfectPartialTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
		faultedPartialTotalDiameter = faultedPartialTotalDiameter /
			(faultedPartialTotalDensity * (grid[Mx] - grid[surfacePos + 1]));
		frankPartialTotalDiameter = frankPartialTotalDiameter /
			(frankPartialTotalDensity * (grid[Mx] - grid[surfacePos + 1]));

		// Set the output precision
		const int outputPrecision = 5;

		// Open the output file
		std::fstream outputFile;
		outputFile.open("Alloy.dat", std::fstream::out | std::fstream::app);
		outputFile << std::setprecision(outputPrecision);

		// Output the data
		outputFile << timestep << " " << time << " " << iTotalDensity << " "
				   << iTotalDiameter << " " << vTotalDensity << " "
				   << vTotalDiameter << " " << voidTotalDensity << " "
				   << voidTotalDiameter << " " << faultedTotalDensity << " "
				   << faultedTotalDiameter << " " << perfectTotalDensity << " "
				   << perfectTotalDiameter << " " << frankTotalDensity << " "
				   << frankTotalDiameter << " " << voidPartialTotalDensity
				   << " " << voidPartialTotalDiameter << " "
				   << faultedPartialTotalDensity << " "
				   << faultedPartialTotalDiameter << " "
				   << perfectPartialTotalDensity << " "
				   << perfectPartialTotalDiameter << " "
				   << frankPartialTotalDensity << " "
				   << frankPartialTotalDiameter << std::endl;

		// Close the output file
		outputFile.close();
	}

	// Restore the PETSC solution array
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorScatter1D")
/**
 * This is a monitoring method that will save 1D plots of the xenon
 * concentration distribution at the middle of the grid.
 */
PetscErrorCode
monitorScatter1D(TS ts, PetscInt timestep, PetscReal time, Vec solution, void*)
{
	perf::ScopedTimer myTimer(scatterTimer);

	// Initial declarations
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi, Mx;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % 200 != 0)
		PetscFunctionReturn(0);

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

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	using NetworkType = core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Region = typename NetworkType::Region;
	auto& network = dynamic_cast<NetworkType&>(solverHandler.getNetwork());
	int networkSize = network.getNumClusters();

	// Get the index of the middle of the grid
	PetscInt ix = Mx / 2;

	// If the middle is on this process
	if (ix >= xs && ix < xs + xm) {
		// Create a DataPoint vector to store the data to give to the data
		// provider for the visualization
		auto myPoints =
			std::make_shared<std::vector<viz::dataprovider::DataPoint>>();

		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[ix];

		for (int i = 0; i < networkSize; i++) {
			// Create a Point with the concentration[i] as the value
			// and add it to myPoints
			auto cluster = network.getCluster(i);
			const Region& clReg = cluster.getRegion();
			for (std::size_t j : makeIntervalRange(clReg[Spec::Xe])) {
				viz::dataprovider::DataPoint aPoint;
				aPoint.value = gridPointSolution[i];
				aPoint.t = time;
				aPoint.x = (double)j;
				myPoints->push_back(aPoint);
			}
		}

		// Get the data provider and give it the points
		scatterPlot1D->getDataProvider()->setDataPoints(myPoints);

		// Change the title of the plot and the name of the data
		std::stringstream title;
		title << "Size Distribution";
		scatterPlot1D->getDataProvider()->setDataName(title.str());
		scatterPlot1D->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		scatterPlot1D->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
					  << "s";
		scatterPlot1D->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "Scatter_TS" << timestep << ".png";
		scatterPlot1D->write(fileName.str());
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorSeries1D")
/**
 * This is a monitoring method that will save 1D plots of many concentrations
 */
PetscErrorCode
monitorSeries1D(TS ts, PetscInt timestep, PetscReal time, Vec solution, void*)
{
	perf::ScopedTimer myTimer(seriesTimer);

	// Initial declarations
	PetscErrorCode ierr;
	const double **solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi;
	double x = 0.0;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % 10 != 0)
		PetscFunctionReturn(0);

	// Get the number of processes
	auto xolotlComm = util::getMPIComm();
	int worldSize;
	MPI_Comm_size(xolotlComm, &worldSize);
	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto& network = solverHandler.getNetwork();
	const int networkSize = network.getNumClusters();

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// To plot a maximum of 18 clusters of the whole benchmark
	const int loopSize = std::min(18, networkSize);

	if (procId == 0) {
		// Create a DataPoint vector to store the data to give to the data
		// provider for the visualization
		std::vector<std::vector<viz::dataprovider::DataPoint>> myPoints(
			loopSize);

		// Loop on the grid
		for (xi = xs; xi < xs + xm; xi++) {
			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[xi];

			for (int i = 0; i < loopSize; i++) {
				// Create a DataPoint with the concentration[i] as the value
				// and add it to myPoints
				viz::dataprovider::DataPoint aPoint;
				aPoint.value = gridPointSolution[i];
				aPoint.t = time;
				aPoint.x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];
				myPoints[i].push_back(aPoint);
			}
		}

		// Loop on the other processes
		for (int i = 1; i < worldSize; i++) {
			// Get the size of the local grid of that process
			int localSize = 0;
			MPI_Recv(
				&localSize, 1, MPI_INT, i, 20, xolotlComm, MPI_STATUS_IGNORE);

			// Loop on their grid
			for (int k = 0; k < localSize; k++) {
				// Get the position
				MPI_Recv(
					&x, 1, MPI_DOUBLE, i, 21, xolotlComm, MPI_STATUS_IGNORE);

				for (int j = 0; j < loopSize; j++) {
					// and the concentrations
					double conc = 0.0;
					MPI_Recv(&conc, 1, MPI_DOUBLE, i, 22, xolotlComm,
						MPI_STATUS_IGNORE);

					// Create a Point with the concentration[i] as the value
					// and add it to myPoints
					viz::dataprovider::DataPoint aPoint;
					aPoint.value = conc;
					aPoint.t = time;
					aPoint.x = x;
					myPoints[j].push_back(aPoint);
				}
			}
		}

		for (int i = 0; i < loopSize; i++) {
			// Get the data provider and give it the points
			auto thePoints =
				std::make_shared<std::vector<viz::dataprovider::DataPoint>>(
					myPoints[i]);
			seriesPlot1D->getDataProvider(i)->setDataPoints(thePoints);
			// TODO: get the name or comp of the cluster
			seriesPlot1D->getDataProvider(i)->setDataName(std::to_string(i));
		}

		// Change the title of the plot
		std::stringstream title;
		title << "Concentrations";
		seriesPlot1D->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		seriesPlot1D->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
					  << "s";
		seriesPlot1D->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "log_series_TS" << timestep << ".ppm";
		seriesPlot1D->write(fileName.str());
	}

	else {
		// Send the value of the local grid size to the master process
		MPI_Send(&xm, 1, MPI_DOUBLE, 0, 20, xolotlComm);

		// Loop on the grid
		for (xi = xs; xi < xs + xm; xi++) {
			// Dump x
			x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];

			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[xi];

			// Send the value of the local position to the master process
			MPI_Send(&x, 1, MPI_DOUBLE, 0, 21, xolotlComm);

			for (int i = 0; i < loopSize; i++) {
				// Send the value of the concentrations to the master process
				MPI_Send(
					&gridPointSolution[i], 1, MPI_DOUBLE, 0, 22, xolotlComm);
			}
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "eventFunction1D")
/**
 * This is a method that checks if the surface should move or bursting happen
 */
PetscErrorCode
eventFunction1D(TS ts, PetscReal time, Vec solution, PetscScalar* fvalue, void*)
{
	perf::ScopedTimer myTimer(eventFuncTimer);

	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi, Mx;
	depthPositions1D.clear();
	fvalue[0] = 1.0, fvalue[1] = 1.0, fvalue[2] = 1.0;

	PetscFunctionBeginUser;

	PetscInt TSNumber = -1;
	ierr = TSGetStepNumber(ts, &TSNumber);

	// Skip if it is the same TS as before
	if (TSNumber == previousTSNumber1D)
		PetscFunctionReturn(0);

	// Set the previous TS number
	previousTSNumber1D = TSNumber;

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

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	CHKERRQ(ierr);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();
	xi = surfacePos + solverHandler.getLeftOffset();

	// Get the network
	auto& network = solverHandler.getNetwork();

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// Get the flux handler to know the flux amplitude.
	auto fluxHandler = solverHandler.getFluxHandler();
	double heliumFluxAmplitude = fluxHandler->getFluxAmplitude();

	// Get the delta time from the previous timestep to this timestep
	double dt = time - solverHandler.getPreviousTime();

	// Work of the moving surface first
	if (solverHandler.moveSurface()) {
		// Write the initial surface position
		if (procId == 0 && util::equal(time, 0.0)) {
			std::ofstream outputFile;
			outputFile.open("surface.txt", std::ios::app);
			outputFile << time << " " << grid[surfacePos + 1] - grid[1]
					   << std::endl;
			outputFile.close();
		}

		// Value to know on which processor is the location of the surface,
		// for MPI usage
		int surfaceProc = 0;

		// if xi is on this process
		if (xi >= xs && xi < xs + xm) {
			// Get the concentrations at xi = surfacePos + 1
			gridPointSolution = solutionArray[xi];

			// Compute the total density of intersitials that escaped from the
			// surface since last timestep using the stored flux
			nInterEvent1D += previousIEventFlux1D * dt;

			// Remove the sputtering yield since last timestep
			nInterEvent1D -= sputteringYield1D * heliumFluxAmplitude * dt;

			// Initialize the value for the flux
			double newFlux = 0.0;

			// Factor for finite difference
			double hxLeft = 0.0, hxRight = 0.0;
			if (xi - 1 >= 0 && xi < Mx) {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			}
			else if (xi - 1 < 0) {
				hxLeft = grid[xi + 1] - grid[xi];
				hxRight = (grid[xi + 2] - grid[xi]) / 2.0;
			}
			else {
				hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
				hxRight = grid[xi + 1] - grid[xi];
			}
			double factor = 2.0 / (hxLeft + hxRight);

			// Consider each interstitial cluster.
			for (int i = 0; i < iClusterIds1D.size(); i++) {
				auto currId = iClusterIds1D[i];
				// Get the cluster
				auto cluster = network.getClusterCommon(currId);
				// Get its concentration
				double conc = gridPointSolution[currId];
				// Get its size and diffusion coefficient
				int size = i + 1;
				double coef = cluster.getDiffusionCoefficient(xi - xs);
				// Compute the flux going to the left
				newFlux += (double)size * factor * coef * conc;
			}

			// Update the previous flux
			previousIEventFlux1D = newFlux;

			// Set the surface processor
			surfaceProc = procId;
		}

		// Get which processor will send the information
		int surfaceId = 0;
		MPI_Allreduce(
			&surfaceProc, &surfaceId, 1, MPI_INT, MPI_SUM, xolotlComm);

		// Send the information about nInterEvent1D and previousFlux1D
		// to the other processes
		MPI_Bcast(&nInterEvent1D, 1, MPI_DOUBLE, surfaceId, xolotlComm);
		MPI_Bcast(&previousIEventFlux1D, 1, MPI_DOUBLE, surfaceId, xolotlComm);

		// Now that all the processes have the same value of nInterstitials,
		// compare it to the threshold to now if we should move the surface

		// Get the initial vacancy concentration
		double initialVConc = solverHandler.getInitialVConc();

		// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
		double threshold = (62.8 - initialVConc) * (grid[xi] - grid[xi - 1]);
		if (nInterEvent1D > threshold) {
			// The surface is moving
			fvalue[0] = 0;
		}

		// Moving the surface back
		else if (nInterEvent1D < -threshold / 10.0) {
			// The surface is moving
			fvalue[1] = 0;
		}
	}

	// Now work on the bubble bursting
	if (solverHandler.burstBubbles()) {
		using NetworkType = core::network::PSIReactionNetwork<
			core::network::PSIFullSpeciesList>;
		using Spec = typename NetworkType::Species;
		auto psiNetwork = dynamic_cast<NetworkType*>(&network);
		auto dof = network.getDOF();

		// Compute the prefactor for the probability (arbitrary)
		double prefactor =
			heliumFluxAmplitude * dt * solverHandler.getBurstingFactor();

		// The depth parameter to know where the bursting should happen
		double depthParam = solverHandler.getTauBursting(); // nm
		// The number of He per V in a bubble
		double heVRatio = solverHandler.getHeVRatio();

		// For now we are not bursting
		bool burst = false;

		// Loop on the full grid of interest
		for (xi = surfacePos + solverHandler.getLeftOffset();
			 xi < Mx - solverHandler.getRightOffset(); xi++) {
			// If this is the locally owned part of the grid
			if (xi >= xs && xi < xs + xm) {
				// Get the distance from the surface
				double distance =
					(grid[xi] + grid[xi + 1]) / 2.0 - grid[surfacePos + 1];

				// Get the pointer to the beginning of the solution data for
				// this grid point
				gridPointSolution = solutionArray[xi];

				using HostUnmanaged = Kokkos::View<double*, Kokkos::HostSpace,
					Kokkos::MemoryUnmanaged>;
				auto hConcs = HostUnmanaged(gridPointSolution, dof);
				auto dConcs = Kokkos::View<double*>("Concentrations", dof);
				deep_copy(dConcs, hConcs);

				// Compute the helium density at this grid point
				double heDensity =
					psiNetwork->getTotalAtomConcentration(dConcs, Spec::He, 1);

				// Compute the radius of the bubble from the number of helium
				double nV = heDensity * (grid[xi + 1] - grid[xi]) / heVRatio;
				double latticeParam = network.getLatticeParameter();
				double tlcCubed = latticeParam * latticeParam * latticeParam;
				double radius = (sqrt(3.0) / 4) * latticeParam +
					cbrt((3.0 * tlcCubed * nV) / (8.0 * core::pi)) -
					cbrt((3.0 * tlcCubed) / (8.0 * core::pi));

				// Add randomness
				double prob = prefactor *
					(1.0 - (distance - radius) / distance) *
					std::min(1.0,
						exp(-(distance - depthParam) / (depthParam * 2.0)));
				double test = solverHandler.getRNG().GetRandomDouble();

				// If the bubble is too big or the probability is high enough
				if (prob > test || radius > distance) {
					burst = true;
					depthPositions1D.push_back(xi);
				}
			}
		}

		// If at least one grid point is bursting
		int localFlag = 1;
		if (burst) {
			// The event is happening
			localFlag = 0;
		}
		// All the processes should call post event
		int flag = -1;
		MPI_Allreduce(&localFlag, &flag, 1, MPI_INT, MPI_MIN, xolotlComm);
		fvalue[2] = flag;
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "postEventFunction1D")
/**
 * This is a method that moves the surface or burst bubbles
 */
PetscErrorCode
postEventFunction1D(TS ts, PetscInt nevents, PetscInt eventList[],
	PetscReal time, Vec solution, PetscBool, void*)
{
	perf::ScopedTimer myTimer(postEventFuncTimer);

	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	PetscInt xs, xm, xi;

	PetscFunctionBeginUser;

	// Check if the surface has moved or a bubble burst
	if (nevents == 0) {
		PetscFunctionReturn(0);
	}

	// Check if both events happened
	if (nevents == 3)
		throw std::string("\nxolotlSolver::Monitor1D: This is not supposed to "
						  "happen, the surface cannot "
						  "move in both directions at the same time!!");

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

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	CHKERRQ(ierr);

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Get the network
	auto& network = solverHandler.getNetwork();
	int dof = network.getDOF();

	// Get the physical grid
	auto grid = solverHandler.getXGrid();

	// Get the flux handler to know the flux amplitude.
	auto fluxHandler = solverHandler.getFluxHandler();
	double heliumFluxAmplitude = fluxHandler->getFluxAmplitude();

	// Get the delta time from the previous timestep to this timestep
	double previousTime = solverHandler.getPreviousTime();
	double dt = time - previousTime;

	// Take care of bursting
	using NetworkType =
		core::network::PSIReactionNetwork<core::network::PSIFullSpeciesList>;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);
	double localNHe = 0.0, localND = 0.0, localNT = 0.0;

	// Loop on each bursting depth
	for (int i = 0; i < depthPositions1D.size(); i++) {
		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[depthPositions1D[i]];

		// Get the distance from the surface
		xi = depthPositions1D[i];
		double distance =
			(grid[xi] + grid[xi + 1]) / 2.0 - grid[surfacePos + 1];
		double hxLeft = 0.0;
		if (xi - 1 < 0) {
			hxLeft = grid[xi + 1] - grid[xi];
		}
		else {
			hxLeft = (grid[xi + 1] - grid[xi - 1]) / 2.0;
		}

		// Write the bursting information
		std::ofstream outputFile;
		outputFile.open("bursting.txt", std::ios::app);
		outputFile << time << " " << distance << std::endl;
		outputFile.close();

		// Pinhole case
		// Loop on every cluster
		for (unsigned int i = 0; i < network.getNumClusters(); i++) {
			const auto& clReg =
				psiNetwork->getCluster(i, plsm::onHost).getRegion();
			// Non-grouped clusters
			if (clReg.isSimplex()) {
				// Get the composition
				Composition comp = clReg.getOrigin();
				// Pure He, D, or T case
				if (comp.isOnAxis(Spec::He)) {
					// Compute the number of atoms released
					localNHe +=
						gridPointSolution[i] * (double)comp[Spec::He] * hxLeft;
					// Reset concentration
					gridPointSolution[i] = 0.0;
				}
				else if (comp.isOnAxis(Spec::D)) {
					// Compute the number of atoms released
					localND +=
						gridPointSolution[i] * (double)comp[Spec::D] * hxLeft;
					// Reset concentration
					gridPointSolution[i] = 0.0;
				}
				else if (comp.isOnAxis(Spec::T)) {
					// Compute the number of atoms released
					localNT +=
						gridPointSolution[i] * (double)comp[Spec::T] * hxLeft;
					// Reset concentration
					gridPointSolution[i] = 0.0;
				}
				// Mixed cluster case
				else if (!comp.isOnAxis(Spec::V) && !comp.isOnAxis(Spec::I)) {
					// Compute the number of atoms released
					localNHe +=
						gridPointSolution[i] * (double)comp[Spec::He] * hxLeft;
					localND +=
						gridPointSolution[i] * (double)comp[Spec::D] * hxLeft;
					localNT +=
						gridPointSolution[i] * (double)comp[Spec::T] * hxLeft;
					// Transfer concentration to V of the same size
					Composition vComp = Composition::zero();
					vComp[Spec::V] = comp[Spec::V];
					auto vCluster =
						psiNetwork->findCluster(vComp, plsm::onHost);
					gridPointSolution[vCluster.getId()] += gridPointSolution[i];
					gridPointSolution[i] = 0.0;
				}
			}
			// Grouped clusters
			else {
				// Compute the number of atoms released
				double concFactor = clReg.volume() / clReg[Spec::He].length();
				for (auto j : makeIntervalRange(clReg[Spec::He])) {
					localNHe +=
						gridPointSolution[i] * (double)j * concFactor * hxLeft;
				}
				concFactor = clReg.volume() / clReg[Spec::D].length();
				for (auto j : makeIntervalRange(clReg[Spec::D])) {
					localND +=
						gridPointSolution[i] * (double)j * concFactor * hxLeft;
				}
				concFactor = clReg.volume() / clReg[Spec::T].length();
				for (auto j : makeIntervalRange(clReg[Spec::T])) {
					localND +=
						gridPointSolution[i] * (double)j * concFactor * hxLeft;
				}

				// Get the factor
				concFactor = clReg.volume() / clReg[Spec::V].length();
				// Loop on the Vs
				for (auto j : makeIntervalRange(clReg[Spec::V])) {
					// Transfer concentration to V of the same size
					Composition vComp = Composition::zero();
					vComp[Spec::V] = j;
					auto vCluster =
						psiNetwork->findCluster(vComp, plsm::onHost);
					// TODO: refine formula with V moment
					gridPointSolution[vCluster.getId()] +=
						gridPointSolution[i] * concFactor;
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

	// Add up the local quantities
	double nHe = 0.0, nD = 0.0, nT = 0.0;
	MPI_Allreduce(&localNHe, &nHe, 1, MPI_DOUBLE, MPI_SUM, xolotlComm);
	nHeliumBurst1D += nHe;
	MPI_Allreduce(&localND, &nD, 1, MPI_DOUBLE, MPI_SUM, xolotlComm);
	nDeuteriumBurst1D += nD;
	MPI_Allreduce(&localNT, &nT, 1, MPI_DOUBLE, MPI_SUM, xolotlComm);
	nTritiumBurst1D += nT;

	// Now takes care of moving surface
	bool moving = false;
	bool movingUp = false;
	for (int i = 0; i < nevents; i++) {
		if (eventList[i] < 2)
			moving = true;
		if (eventList[i] == 0)
			movingUp = true;
	}

	// Skip if nothing is moving
	if (!moving) {
		// Restore the solutionArray
		ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
		CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}

	// Set the surface position
	xi = surfacePos + solverHandler.getLeftOffset();

	// Get the initial vacancy concentration
	double initialVConc = solverHandler.getInitialVConc();

	// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
	double threshold = (62.8 - initialVConc) * (grid[xi] - grid[xi - 1]);

	if (movingUp) {
		int nGridPoints = 0;
		// Move the surface up until it is smaller than the next threshold
		while (nInterEvent1D > threshold) {
			// Move the surface higher
			surfacePos--;
			xi = surfacePos + solverHandler.getLeftOffset();
			nGridPoints++;
			// Update the number of interstitials
			nInterEvent1D -= threshold;
			// Update the thresold
			threshold = (62.8 - initialVConc) * (grid[xi] - grid[xi - 1]);
		}

		// Throw an exception if the position is negative
		if (surfacePos < 0) {
			PetscBool flagCheck;
			ierr =
				PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCheck);
			CHKERRQ(ierr);
			if (flagCheck) {
				// Write the convergence reason
				std::ofstream outputFile;
				outputFile.open("solverStatus.txt");
				outputFile << "overgrid" << std::endl;
				outputFile.close();
			}
			throw std::string("\nxolotlSolver::Monitor1D: The surface is "
							  "trying to go outside of the grid!!");
		}

		// Printing information about the extension of the material
		if (procId == 0) {
			std::cout << "Adding " << nGridPoints
					  << " points to the grid at time: " << time << " s."
					  << std::endl;
		}

		// Set it in the solver
		solverHandler.setSurfacePosition(surfacePos);

		// Initialize the vacancy concentration and the temperature on the new
		// grid points Get the single vacancy ID
		auto singleVacancyCluster = network.getSingleVacancy();
		auto vacancyIndex = core::network::IReactionNetwork::invalidIndex();
		if (singleVacancyCluster.getId() !=
			core::network::IReactionNetwork::invalidIndex())
			vacancyIndex = singleVacancyCluster.getId();
		// Get the surface temperature
		double temp = 0.0;
		if (xi >= xs && xi < xs + xm) {
			temp = solutionArray[xi][dof];
		}
		double surfTemp = 0.0;
		MPI_Allreduce(&temp, &surfTemp, 1, MPI_DOUBLE, MPI_SUM, xolotlComm);

		// Loop on the new grid points
		while (nGridPoints >= 0) {
			// Position of the newly created grid point
			xi = surfacePos + nGridPoints;

			// If xi is on this process
			if (xi >= xs && xi < xs + xm) {
				// Get the concentrations
				gridPointSolution = solutionArray[xi];

				// Set the new surface temperature
				gridPointSolution[dof] = surfTemp;

				if (vacancyIndex !=
						core::network::IReactionNetwork::invalidIndex() &&
					nGridPoints > 0) {
					// Initialize the vacancy concentration
					gridPointSolution[vacancyIndex] = initialVConc;
				}
			}

			// Decrease the number of grid points
			--nGridPoints;
		}
	}

	// Moving the surface back
	else {
		// Move it back as long as the number of interstitials in negative
		while (nInterEvent1D < 0.0) {
			// Compute the threshold to a deeper grid point
			threshold = (62.8 - initialVConc) * (grid[xi + 1] - grid[xi]);
			// Set all the concentrations to 0.0 at xi = surfacePos + 1
			// if xi is on this process
			if (xi >= xs && xi < xs + xm) {
				// Get the concentrations at xi = surfacePos + 1
				gridPointSolution = solutionArray[xi];
				// Loop on DOF
				for (int i = 0; i < dof; i++) {
					gridPointSolution[i] = 0.0;
				}
			}

			// Move the surface deeper
			surfacePos++;
			xi = surfacePos + solverHandler.getLeftOffset();
			// Update the number of interstitials
			nInterEvent1D += threshold;
		}

		// Printing information about the extension of the material
		if (procId == 0) {
			std::cout << "Removing grid points to the grid at time: " << time
					  << " s." << std::endl;
		}

		// Set it in the solver
		solverHandler.setSurfacePosition(surfacePos);
	}

	// Set the new surface location in the surface advection handler
	auto advecHandler = solverHandler.getAdvectionHandler();
	advecHandler->setLocation(grid[surfacePos + 1] - grid[1]);

	// Set the new surface in the temperature handler
	auto tempHandler = solverHandler.getTemperatureHandler();
	tempHandler->updateSurfacePosition(surfacePos);

	// Get the flux handler to reinitialize it
	fluxHandler->initializeFluxHandler(
		solverHandler.getNetwork(), surfacePos, grid);

	// Get the modified trap-mutation handler to reinitialize it
	auto mutationHandler = solverHandler.getMutationHandler();
	auto advecHandlers = solverHandler.getAdvectionHandlers();
	mutationHandler->initializeIndex1D(
		surfacePos, solverHandler.getNetwork(), advecHandlers, grid, xm, xs);

	// Write the updated surface position
	if (procId == 0) {
		std::ofstream outputFile;
		outputFile.open("surface.txt", std::ios::app);
		outputFile << time << " " << grid[surfacePos + 1] - grid[1]
				   << std::endl;
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
setupPetsc1DMonitor(TS ts)
{
	PetscErrorCode ierr;

	auto handlerRegistry = perf::getHandlerRegistry();

	// Initialize the timers, including the one for this function.
	initTimer = handlerRegistry->getTimer("monitor1D:init");
	perf::ScopedTimer myTimer(initTimer);
	checkNegativeTimer = handlerRegistry->getTimer("monitor1D:checkNeg");
	tridynTimer = handlerRegistry->getTimer("monitor1D:tridyn");
	startStopTimer = handlerRegistry->getTimer("monitor1D:startStop");
	heRetentionTimer = handlerRegistry->getTimer("monitor1D:heRet");
	xeRetentionTimer = handlerRegistry->getTimer("monitor1D:xeRet");
	scatterTimer = handlerRegistry->getTimer("monitor1D:scatter");
	seriesTimer = handlerRegistry->getTimer("monitor1D:series");
	eventFuncTimer = handlerRegistry->getTimer("monitor1D:event");
	postEventFuncTimer = handlerRegistry->getTimer("monitor1D:postEvent");

	// Get the process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get xolotlViz handler registry
	auto vizHandlerRegistry = viz::VizHandlerRegistry::get();

	// Flags to launch the monitors or not
	PetscBool flagNeg, flagCollapse, flag2DPlot, flag1DPlot, flagSeries,
		flagPerf, flagHeRetention, flagStatus, flagXeRetention, flagTRIDYN,
		flagAlloy, flagTemp, flagLargest;

	// Check the option -check_negative
	ierr = PetscOptionsHasName(NULL, NULL, "-check_negative", &flagNeg);
	checkPetscError(ierr,
		"setupPetsc1DMonitor: PetscOptionsHasName (-check_negative) failed.");

	// Check the option -check_collapse
	ierr = PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCollapse);
	checkPetscError(ierr,
		"setupPetsc1DMonitor: PetscOptionsHasName (-check_collapse) failed.");

	// Check the option -plot_perf
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_perf", &flagPerf);
	checkPetscError(
		ierr, "setupPetsc1DMonitor: PetscOptionsHasName (-plot_perf) failed.");

	// Check the option -plot_series
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_series", &flagSeries);
	checkPetscError(ierr,
		"setupPetsc1DMonitor: PetscOptionsHasName (-plot_series) failed.");

	// Check the option -plot_1d
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_1d", &flag1DPlot);
	checkPetscError(
		ierr, "setupPetsc1DMonitor: PetscOptionsHasName (-plot_1d) failed.");

	// Check the option -helium_retention
	ierr =
		PetscOptionsHasName(NULL, NULL, "-helium_retention", &flagHeRetention);
	checkPetscError(ierr,
		"setupPetsc1DMonitor: PetscOptionsHasName (-helium_retention) failed.");

	// Check the option -xenon_retention
	ierr =
		PetscOptionsHasName(NULL, NULL, "-xenon_retention", &flagXeRetention);
	checkPetscError(ierr,
		"setupPetsc1DMonitor: PetscOptionsHasName (-xenon_retention) failed.");

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, NULL, "-start_stop", &flagStatus);
	checkPetscError(
		ierr, "setupPetsc1DMonitor: PetscOptionsHasName (-start_stop) failed.");

	// Check the option -tridyn
	ierr = PetscOptionsHasName(NULL, NULL, "-tridyn", &flagTRIDYN);
	checkPetscError(
		ierr, "setupPetsc1DMonitor: PetscOptionsHasName (-tridyn) failed.");

	// Check the option -alloy
	ierr = PetscOptionsHasName(NULL, NULL, "-alloy", &flagAlloy);
	checkPetscError(
		ierr, "setupPetsc1DMonitor: PetscOptionsHasName (-alloy) failed.");

	// Check the option -temp_profile
	ierr = PetscOptionsHasName(NULL, NULL, "-temp_profile", &flagTemp);
	checkPetscError(ierr,
		"setupPetsc1DMonitor: PetscOptionsHasName (-temp_profile) failed.");

	// Check the option -largest_conc
	ierr = PetscOptionsHasName(NULL, NULL, "-largest_conc", &flagLargest);
	checkPetscError(ierr,
		"setupPetsc1DMonitor: PetscOptionsHasName (-largest_conc) failed.");

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto& network = solverHandler.getNetwork();
	const int networkSize = network.getNumClusters();
	// Get the number of species
	auto numSpecies = network.getSpeciesListSize();

	// Create data depending on the boundary conditions
	if (solverHandler.getLeftOffset() == 1) {
		nSurf1D = std::vector<double>(numSpecies, 0.0);
		previousSurfFlux1D = std::vector<double>(numSpecies, 0.0);
	}
	if (solverHandler.getRightOffset() == 1) {
		nBulk1D = std::vector<double>(numSpecies, 0.0);
		previousBulkFlux1D = std::vector<double>(numSpecies, 0.0);
	}

	// Determine if we have an existing restart file,
	// and if so, it it has had timesteps written to it.
	std::unique_ptr<io::XFile> networkFile;
	std::unique_ptr<io::XFile::TimestepGroup> lastTsGroup;
	std::string networkName = solverHandler.getNetworkName();
	bool hasConcentrations = false;
	if (not networkName.empty()) {
		networkFile.reset(new io::XFile(networkName));
		auto concGroup = networkFile->getGroup<io::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
		if (hasConcentrations) {
			lastTsGroup = concGroup->getLastTimestepGroup();
		}
	}

	// Set the post step processing to stop the solver if the time step
	// collapses
	if (flagCollapse) {
		// Find the threshold
		PetscBool flag;
		ierr = PetscOptionsGetReal(
			NULL, NULL, "-check_collapse", &timeStepThreshold, &flag);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsGetReal (-check_collapse) "
			"failed.");
		if (!flag)
			timeStepThreshold = 1.0e-16;

		// Set the post step process that tells the solver when to stop if the
		// time step collapse
		ierr = TSSetPostStep(ts, checkTimeStep);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSSetPostStep (checkTimeStep) failed.");
	}

	// Set the monitor to check the negative concentrations
	if (flagNeg) {
		// Find the stride to know how often we want to check
		PetscBool flag;
		ierr = PetscOptionsGetReal(
			NULL, NULL, "-check_negative", &negThreshold1D, &flag);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsGetReal (-check_negative) "
			"failed.");
		if (!flag)
			negThreshold1D = 1.0e-30;

		// checkNegative1D will be called at each timestep
		ierr = TSMonitorSet(ts, checkNegative1D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (checkNegative1D) failed.");
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr = PetscOptionsGetReal(
			NULL, NULL, "-start_stop", &hdf5Stride1D, &flag);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsGetReal (-start_stop) failed.");
		if (!flag)
			hdf5Stride1D = 1.0;

		// Compute the correct hdf5Previous1D for a restart
		// Get the last time step written in the HDF5 file
		if (hasConcentrations) {
			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			solverHandler.setPreviousTime(previousTime);
			hdf5Previous1D = (int)(previousTime / hdf5Stride1D);
		}

		// Don't do anything if both files have the same name
		if (hdf5OutputName1D != solverHandler.getNetworkName()) {
			PetscInt Mx;
			PetscErrorCode ierr;

			// Get the da from ts
			DM da;
			ierr = TSGetDM(ts, &da);
			checkPetscError(ierr, "setupPetsc1DMonitor: TSGetDM failed.");

			// Get the size of the total grid
			ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE,
				PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
				PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
				PETSC_IGNORE, PETSC_IGNORE);
			checkPetscError(ierr, "setupPetsc1DMonitor: DMDAGetInfo failed.");

			// Get the solver handler
			auto& solverHandler = PetscSolver::getSolverHandler();

			// Get the physical grid
			auto grid = solverHandler.getXGrid();

			// Create and initialize a checkpoint file.
			// We do this in its own scope so that the file
			// is closed when the file object goes out of scope.
			// We want it to close before we (potentially) copy
			// the network from another file using a single-process
			// MPI communicator.
			{
				io::XFile checkpointFile(hdf5OutputName1D, grid, xolotlComm);
			}

			// Copy the network group from the given file (if it has one).
			// We open the files using a single-process MPI communicator
			// because it is faster for a single process to do the
			// copy with HDF5's H5Ocopy implementation than it is
			// when all processes call the copy function.
			// The checkpoint file must be closed before doing this.
			writeNetwork(xolotlComm, solverHandler.getNetworkName(),
				hdf5OutputName1D, network);
		}

		// startStop1D will be called at each timestep
		ierr = TSMonitorSet(ts, startStop1D, NULL, NULL);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (startStop1D) failed.");
	}

	// If the user wants the surface to be able to move or bursting
	if (solverHandler.moveSurface() || solverHandler.burstBubbles()) {
		// Surface
		if (solverHandler.moveSurface()) {
			using NetworkType = core::network::PSIReactionNetwork<
				core::network::PSIFullSpeciesList>;
			auto psiNetwork = dynamic_cast<NetworkType*>(&network);

			// Initialize the composition
			NetworkType::Composition comp = NetworkType::Composition::zero();

			// Loop on interstital clusters
			bool iClusterExists = true;
			int iSize = 1;
			while (iClusterExists) {
				comp[NetworkType::Species::I] = iSize;
				auto cluster = psiNetwork->findCluster(comp, plsm::onHost);
				// Check that the helium cluster is present in the network
				if (cluster.getId() != NetworkType::invalidIndex()) {
					iClusterIds1D.push_back(cluster.getId());
					iSize++;
				}
				else
					iClusterExists = false;
			}

			// Get the interstitial information at the surface if concentrations
			// were stored
			if (hasConcentrations) {
				assert(lastTsGroup);

				// Get the interstitial quantity from the HDF5 file
				nInterEvent1D = lastTsGroup->readData1D("nInterstitial");
				// Get the previous I flux from the HDF5 file
				previousIEventFlux1D = lastTsGroup->readData1D("previousIFlux");
				// Get the previous time from the HDF5 file
				double previousTime = lastTsGroup->readPreviousTime();
				solverHandler.setPreviousTime(previousTime);
			}

			// Get the sputtering yield
			sputteringYield1D = solverHandler.getSputteringYield();

			// Master process
			if (procId == 0) {
				// Clear the file where the surface will be written
				std::ofstream outputFile;
				outputFile.open("surface.txt");
				outputFile.close();
			}
		}

		// Set directions and terminate flags for the surface event
		PetscInt direction[3];
		PetscBool terminate[3];
		direction[0] = 0, direction[1] = 0, direction[2] = 0;
		terminate[0] = PETSC_FALSE, terminate[1] = PETSC_FALSE,
		terminate[2] = PETSC_FALSE;
		// Set the TSEvent
		ierr = TSSetEventHandler(ts, 3, direction, terminate, eventFunction1D,
			postEventFunction1D, NULL);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSSetEventHandler (eventFunction1D) failed.");

		if (solverHandler.burstBubbles() && procId == 0) {
			// Uncomment to clear the file where the bursting info will be
			// written
			std::ofstream outputFile;
			outputFile.open("bursting.txt");
			outputFile.close();
		}
	}

	// Set the monitor to save 1D plot of xenon distribution
	if (flag1DPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a ScatterPlot
			scatterPlot1D = vizHandlerRegistry->getPlot(
				"scatterPlot1D", viz::PlotType::SCATTER);

			scatterPlot1D->setLogScale();

			// Create and set the label provider
			auto labelProvider =
				std::make_shared<viz::LabelProvider>("labelProvider");
			labelProvider->axis1Label = "Xenon Size";
			labelProvider->axis2Label = "Concentration";

			// Give it to the plot
			scatterPlot1D->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider =
				std::make_shared<viz::dataprovider::CvsXDataProvider>(
					"dataProvider");

			// Give it to the plot
			scatterPlot1D->setDataProvider(dataProvider);
		}

		// monitorScatter1D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorScatter1D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (monitorScatter1D) failed.");
	}

	// Set the monitor to save 1D plot of many concentrations
	if (flagSeries) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a ScatterPlot
			seriesPlot1D = vizHandlerRegistry->getPlot(
				"seriesPlot1D", viz::PlotType::SERIES);

			// set the log scale
			//			seriesPlot1D->setLogScale();

			// Create and set the label provider
			auto labelProvider =
				std::make_shared<viz::LabelProvider>("labelProvider");
			labelProvider->axis1Label = "x Position on the Grid";
			labelProvider->axis2Label = "Concentration";

			// Give it to the plot
			seriesPlot1D->setLabelProvider(labelProvider);

			// To plot a maximum of 18 clusters of the whole benchmark
			const int loopSize = std::min(18, networkSize);

			// Create a data provider for each cluster in the network
			for (int i = 0; i < loopSize; i++) {
				// Set the name for Identifiable
				std::stringstream dataProviderName;
				dataProviderName << "dataprovider" << i;
				// Create the data provider
				auto dataProvider =
					std::make_shared<viz::dataprovider::CvsXDataProvider>(
						dataProviderName.str());

				// Give it to the plot
				seriesPlot1D->addDataProvider(dataProvider);
			}
		}

		// monitorSeries1D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSeries1D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (monitorSeries1D) failed.");
	}

	// Set the monitor to save performance plots (has to be in parallel)
	if (flagPerf) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a ScatterPlot
			perfPlot =
				vizHandlerRegistry->getPlot("perfPlot", viz::PlotType::SCATTER);

			// Create and set the label provider
			auto labelProvider =
				std::make_shared<viz::LabelProvider>("labelProvider");
			labelProvider->axis1Label = "Process ID";
			labelProvider->axis2Label = "Solver Time";

			// Give it to the plot
			perfPlot->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider =
				std::make_shared<viz::dataprovider::CvsXDataProvider>(
					"dataProvider");

			// Give it to the plot
			perfPlot->setDataProvider(dataProvider);
		}

		// monitorPerf will be called at each timestep
		ierr = TSMonitorSet(ts, monitorPerf, NULL, NULL);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (monitorPerf) failed.");
	}

	// Set the monitor to compute the helium fluence and the retention
	// for the retention calculation
	if (flagHeRetention) {
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

			// Get the names of the species in the network
			std::vector<std::string> names;
			for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
				names.push_back(network.getSpeciesName(id));
			}

			// If the surface is a free surface
			if (solverHandler.getLeftOffset() == 1) {
				// Loop on the names
				for (auto i = 0; i < names.size(); i++) {
					// Create the n attribute name
					std::ostringstream nName;
					nName << "n" << names[i] << "Surf";
					// Read quantity attribute
					nSurf1D[i] = lastTsGroup->readData1D(nName.str());

					// Create the previous flux attribute name
					std::ostringstream prevFluxName;
					prevFluxName << "previousFlux" << names[i] << "Surf";
					// Read the attribute
					previousSurfFlux1D[i] =
						lastTsGroup->readData1D(prevFluxName.str());
				}
			}

			// If the bottom is a free surface
			if (solverHandler.getRightOffset() == 1) {
				// Loop on the names
				for (auto i = 0; i < names.size(); i++) {
					// Create the n attribute name
					std::ostringstream nName;
					nName << "n" << names[i] << "Bulk";
					// Read quantity attribute
					nBulk1D[i] = lastTsGroup->readData1D(nName.str());

					// Create the previous flux attribute name
					std::ostringstream prevFluxName;
					prevFluxName << "previousFlux" << names[i] << "Bulk";
					// Read the attribute
					previousBulkFlux1D[i] =
						lastTsGroup->readData1D(prevFluxName.str());
				}
			}

			// Bursting
			if (solverHandler.burstBubbles()) {
				// Read about the impurity fluxes in from bursting
				nHeliumBurst1D = lastTsGroup->readData1D("nHeliumBurst");
				nDeuteriumBurst1D = lastTsGroup->readData1D("nDeuteriumBurst");
				nTritiumBurst1D = lastTsGroup->readData1D("nTritiumBurst");
			}
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeFluence, NULL, NULL);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeHeliumRetention1D will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumRetention1D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (computeHeliumRetention1D) "
			"failed.");

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
		// Get the da from ts
		DM da;
		ierr = TSGetDM(ts, &da);
		checkPetscError(ierr, "setupPetsc1DMonitor: TSGetDM failed.");
		// Get the local boundaries
		PetscInt xm;
		ierr = DMDAGetCorners(da, NULL, NULL, NULL, &xm, NULL, NULL);
		checkPetscError(ierr, "setupPetsc1DMonitor: DMDAGetCorners failed.");
		// Create the local vectors on each process
		solverHandler.createLocalNE(xm);

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
			ierr, "setupPetsc1DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeXenonRetention1D will be called at each timestep
		ierr = TSMonitorSet(ts, computeXenonRetention1D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (computeXenonRetention1D) "
			"failed.");

		// Master process
		if (procId == 0) {
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("retentionOut.txt");
			outputFile.close();
		}
	}

	// Set the monitor to output data for TRIDYN
	if (flagTRIDYN) {
		// computeTRIDYN1D will be called at each timestep
		ierr = TSMonitorSet(ts, computeTRIDYN1D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (computeTRIDYN1D) failed.");
	}

	// Set the monitor to output data for Alloy
	if (flagAlloy) {
		if (procId == 0) {
			// Create/open the output files
			std::fstream outputFile;
			outputFile.open("Alloy.dat", std::fstream::out);
			outputFile.close();
		}

		// computeAlloy1D will be called at each timestep
		ierr = TSMonitorSet(ts, computeAlloy1D, NULL, NULL);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (computeAlloy1D) failed.");
	}

	// Set the monitor to compute the temperature profile
	if (flagTemp) {
		if (procId == 0) {
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("tempProf.txt");

			// Get the da from ts
			DM da;
			ierr = TSGetDM(ts, &da);
			checkPetscError(ierr, "setupPetsc1DMonitor: TSGetDM failed.");

			// Get the total size of the grid
			PetscInt Mx;
			ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE,
				PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
				PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
				PETSC_IGNORE, PETSC_IGNORE);
			checkPetscError(ierr, "setupPetsc1DMonitor: DMDAGetInfo failed.");

			// Get the physical grid
			auto grid = solverHandler.getXGrid();
			// Get the position of the surface
			int surfacePos = solverHandler.getSurfacePosition();

			// Loop on the entire grid
			for (int xi = surfacePos + solverHandler.getLeftOffset();
				 xi < Mx - solverHandler.getRightOffset(); xi++) {
				// Set x
				double x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];
				outputFile << x << " ";
			}
			outputFile << std::endl;
			outputFile.close();
		}

		// computeCumulativeHelium1D will be called at each timestep
		ierr = TSMonitorSet(ts, profileTemperature1D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (profileTemperature1D) failed.");
	}

	// Set the monitor to monitor the concentration of the largest cluster
	if (flagLargest) {
		// Look for the largest cluster
		auto& network = solverHandler.getNetwork();
		largestClusterId1D = network.getLargestClusterId();

		// Find the threshold
		PetscBool flag;
		ierr = PetscOptionsGetReal(
			NULL, NULL, "-largest_conc", &largestThreshold1D, &flag);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsGetReal (-largest_conc) failed.");

		// monitorLargest1D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorLargest1D, NULL, NULL);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (monitorLargest1D) failed.");
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(ts, monitorTime, NULL, NULL);
	checkPetscError(
		ierr, "setupPetsc1DMonitor: TSMonitorSet (monitorTime) failed.");

	PetscFunctionReturn(0);
}

} /* end namespace monitor */
} /* end namespace solver */
} /* end namespace xolotl */
