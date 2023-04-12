#include <xolotl/core/network/AlloyReactionNetwork.h>
#include <xolotl/core/network/FeReactionNetwork.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/core/network/ZrReactionNetwork.h>
#include <xolotl/io/XFile.h>
#include <xolotl/solver/PetscSolver.h>
#include <xolotl/solver/monitor/PetscMonitor0D.h>
#include <xolotl/solver/monitor/PetscMonitorFunctions.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/viz/dataprovider/CvsXDataProvider.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
/**
 * This is a monitoring method that will compute average density and diameter
 * of defects.
 */
PetscErrorCode
computeAlphaZr(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<PetscMonitor0D*>(ictx)->computeAlphaZr(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
monitorBubble(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<PetscMonitor0D*>(ictx)->monitorBubble(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

void
PetscMonitor0D::setup(int loop)
{
	PetscErrorCode ierr;

	_loopNumber = loop;

	// Get xolotlViz handler registry
	auto vizHandlerRegistry = _solverHandler->getVizHandler();

	// Flags to launch the monitors or not
	PetscBool flagCheck, flag1DPlot, flagBubble, flagStatus, flagAlloy,
		flagXeRetention, flagHeRetention, flagLargest, flagZr;

	// Check the option -check_collapse
	ierr = PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCheck);
	checkPetscError(ierr,
		"setupPetsc0DMonitor: PetscOptionsHasName (-check_collapse) failed.");

	// Check the option -plot_1d
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_1d", &flag1DPlot);
	checkPetscError(
		ierr, "setupPetsc0DMonitor: PetscOptionsHasName (-plot_1d) failed.");

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, NULL, "-start_stop", &flagStatus);
	checkPetscError(
		ierr, "setupPetsc0DMonitor: PetscOptionsHasName (-start_stop) failed.");

	// Check the option -bubble
	ierr = PetscOptionsHasName(NULL, NULL, "-bubble", &flagBubble);
	checkPetscError(
		ierr, "setupPetsc0DMonitor: PetscOptionsHasName (-bubble) failed.");

	// Check the option -alloy
	ierr = PetscOptionsHasName(NULL, NULL, "-alloy", &flagAlloy);
	checkPetscError(
		ierr, "setupPetsc0DMonitor: PetscOptionsHasName (-alloy) failed.");

	// Check the option -alpha_zr
	ierr = PetscOptionsHasName(NULL, NULL, "-alpha_zr", &flagZr);
	checkPetscError(
		ierr, "setupPetsc0DMonitor: PetscOptionsHasName (-alpha_zr) failed.");

	// Check the option -xenon_retention
	ierr =
		PetscOptionsHasName(NULL, NULL, "-xenon_retention", &flagXeRetention);
	checkPetscError(ierr,
		"setupPetsc0DMonitor: PetscOptionsHasName (-xenon_retention) failed.");

	// Check the option -helium_retention
	ierr =
		PetscOptionsHasName(NULL, NULL, "-helium_retention", &flagHeRetention);
	checkPetscError(ierr,
		"setupPetsc0DMonitor: PetscOptionsHasName (-helium_retention) failed.");

	// Check the option -largest_conc
	ierr = PetscOptionsHasName(NULL, NULL, "-largest_conc", &flagLargest);
	checkPetscError(ierr,
		"setupPetsc0DMonitor: PetscOptionsHasName (-largest_conc) failed.");

	// Determine if we have an existing restart file,
	// and if so, it it has had timesteps written to it.
	std::unique_ptr<io::XFile> networkFile;
	std::unique_ptr<io::XFile::TimestepGroup> lastTsGroup;
	std::string networkName = _solverHandler->getNetworkName();
	bool hasConcentrations = false;
	if (not networkName.empty()) {
		networkFile = std::make_unique<io::XFile>(networkName);
		auto concGroup = networkFile->getGroup<io::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
		if (hasConcentrations) {
			lastTsGroup = concGroup->getLastTimestepGroup();
		}
	}

	// Set the post step processing to stop the solver if the time step
	// collapses
	if (flagCheck) {
		// Find the threshold
		PetscBool flag;
		ierr = PetscOptionsGetReal(
			NULL, NULL, "-check_collapse", &timeStepThreshold, &flag);
		checkPetscError(ierr,
			"setupPetsc0DMonitor: PetscOptionsGetInt (-check_collapse) "
			"failed.");
		if (!flag)
			timeStepThreshold = 1.0e-16;

		// Set the post step process that tells the solver when to stop if the
		// time step collapse
		ierr = TSSetPostStep(_ts, checkTimeStep);
		checkPetscError(
			ierr, "setupPetsc0DMonitor: TSSetPostStep (checkTimeStep) failed.");
	}

	// Set the monitor to save 1D plot of xenon distribution
	if (flag1DPlot) {
		// Create a ScatterPlot
		_scatterPlot = vizHandlerRegistry->getPlot(viz::PlotType::SCATTER);

		//		_scatterPlot->setLogScale();

		// Create and set the label provider
		auto labelProvider = std::make_shared<viz::LabelProvider>();
		labelProvider->axis1Label = "Xenon Size";
		labelProvider->axis2Label = "Concentration";

		// Give it to the plot
		_scatterPlot->setLabelProvider(labelProvider);

		// Create the data provider
		auto dataProvider =
			std::make_shared<viz::dataprovider::CvsXDataProvider>();

		// Give it to the plot
		_scatterPlot->setDataProvider(dataProvider);

		// monitorScatter will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorScatter, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc0DMonitor: TSMonitorSet (monitorScatter) failed.");
	}

	// Set the monitor to save text file of the mean concentration of bubbles
	if (flagBubble) {
		// monitorBubble0D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorBubble, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc0DMonitor: TSMonitorSet (monitorBubble) failed.");
	}

	// Set the monitor to output data for Alloy
	if (flagAlloy) {
		auto& network = _solverHandler->getNetwork();
		auto numSpecies = network.getSpeciesListSize();
		// Create/open the output files
		std::fstream outputFile;
		outputFile.open("Alloy.dat", std::fstream::out);
		outputFile << "#time_step time ";
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			auto speciesName = network.getSpeciesName(id);
			outputFile << speciesName << "_density " << speciesName
					   << "_diameter " << speciesName << "_partial_density "
					   << speciesName << "_partial_diameter ";
		}
		outputFile << std::endl;
		outputFile.close();

		// computeAlloy0D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeAlloy, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc0DMonitor: TSMonitorSet (computeAlloy) failed.");
	}
	// Set the monitor to output data for AlphaZr
	if (flagZr) {
		auto& network = _solverHandler->getNetwork();
		auto numSpecies = network.getSpeciesListSize();

		// Create/open the output files
		std::fstream outputFile;
		outputFile.open("AlphaZr.dat", std::fstream::out);
		outputFile << "#time_step time ";
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			auto speciesName = network.getSpeciesName(id);
			outputFile << speciesName << "_density " << speciesName << "_atom "
					   << speciesName << "_diameter " << speciesName
					   << "_partial_density " << speciesName << "_partial_atom "
					   << speciesName << "_partial_diameter ";
		}

		outputFile << std::endl;
		outputFile.close();

		// computeAlphaZr will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeAlphaZr, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc0DMonitor: TSMonitorSet (computeAlphaZr) failed.");
	}

	// Set the monitor to compute the xenon content
	if (flagXeRetention) {
		// Get the previous time if concentrations were stored and initialize
		// the fluence
		if (hasConcentrations) {
			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			_solverHandler->setPreviousTime(previousTime);
			// Initialize the fluence
			auto fluxHandler = _solverHandler->getFluxHandler();
			// Increment the fluence with the value at this current timestep
			fluxHandler->computeFluence(previousTime);
		}

		// computeXenonRetention0D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeXenonRetention, this, nullptr);
		checkPetscError(ierr,
			"setupPetsc0DMonitor: TSMonitorSet (computeXenonRetention) "
			"failed.");

		// Uncomment to clear the file where the retention will be written
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt");
		outputFile << "#time Xenon_conc radius partial_radius "
					  "partial_bubble_conc partial_size"
				   << std::endl;
		outputFile.close();
	}

	// Set the monitor to compute the helium content
	if (flagHeRetention) {
		// Get the previous time if concentrations were stored and initialize
		// the fluence
		if (hasConcentrations) {
			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			_solverHandler->setPreviousTime(previousTime);
			// Initialize the fluence
			auto fluxHandler = _solverHandler->getFluxHandler();
			// Increment the fluence with the value at this current timestep
			fluxHandler->computeFluence(previousTime);
		}

		// computeHeliumRetention0D will be called at each timestep
		ierr =
			TSMonitorSet(_ts, monitor::computeHeliumRetention, this, nullptr);
		checkPetscError(ierr,
			"setupPetsc0DMonitor: TSMonitorSet (computeHeliumRetention) "
			"failed.");

		// Uncomment to clear the file where the retention will be written
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt");
		outputFile << "#time Helium_conc C_b av_He av_V" << std::endl;
		outputFile.close();
	}

	// Set the monitor to monitor the concentration of the largest cluster
	if (flagLargest) {
		// Look for the largest cluster
		auto& network = _solverHandler->getNetwork();
		_largestClusterId = network.getLargestClusterId();

		// Find the threshold
		PetscBool flag;
		ierr = PetscOptionsGetReal(
			NULL, NULL, "-largest_conc", &_largestThreshold, &flag);
		checkPetscError(ierr,
			"setupPetsc0DMonitor: PetscOptionsGetReal (-largest_conc) failed.");

		// monitorLargest1D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorLargest, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc0DMonitor: TSMonitorSet (monitorLargest) failed.");
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr =
			PetscOptionsGetReal(NULL, NULL, "-start_stop", &_hdf5Stride, &flag);
		checkPetscError(ierr,
			"setupPetsc0DMonitor: PetscOptionsGetInt (-start_stop) failed.");
		if (!flag)
			_hdf5Stride = 1.0;

		// Compute the correct _hdf5Previous for a restart
		// Get the last time step written in the HDF5 file
		if (hasConcentrations) {
			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			_solverHandler->setPreviousTime(previousTime);
			_hdf5Previous = (PetscInt)(previousTime / _hdf5Stride);
		}

		// Don't do anything if both files have the same name
		if (_hdf5OutputName != _solverHandler->getNetworkName()) {
			PetscErrorCode ierr;

			// Get the network
			auto& network = _solverHandler->getNetwork();

			// Get the physical grid (which is empty)
			auto grid = _solverHandler->getXGrid();

			// Get the MPI communicator
			auto xolotlComm = util::getMPIComm();

			// Create and initialize a checkpoint file.
			// We do this in its own scope so that the file
			// is closed when the file object goes out of scope.
			// We want it to close before we (potentially) copy
			// the network from another file using a single-process
			// MPI communicator.
			{
				io::XFile checkpointFile(_hdf5OutputName, 1, xolotlComm);
			}

			// Copy the network group from the given file (if it has one).
			// We open the files using a single-process MPI communicator
			// because it is faster for a single process to do the
			// copy with HDF5's H5Ocopy implementation than it is
			// when all processes call the copy function.
			// The checkpoint file must be closed before doing this.
			writeNetwork(
				xolotlComm, _hdf5OutputName, _solverHandler->getNetworkName());
		}

		// startStop0D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::startStop, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc0DMonitor: TSMonitorSet (startStop) failed.");
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(_ts, monitor::monitorTime, this, nullptr);
	checkPetscError(
		ierr, "setupPetsc0DMonitor: TSMonitorSet (monitorTime) failed.");
}

PetscErrorCode
PetscMonitor0D::monitorLargest(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;

	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];
	// Check the concentration
	if (gridPointSolution[_largestClusterId] > _largestThreshold) {
		ierr = TSSetConvergedReason(ts, TS_CONVERGED_USER);
		CHKERRQ(ierr);
		// Send an error
		throw std::runtime_error(
			"\nxolotlSolver::Monitor0D: The largest cluster "
			"concentration is too high!!");
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::startStop(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declaration
	PetscErrorCode ierr;
	const double **solutionArray, *gridPointSolution;

	PetscFunctionBeginUser;

	// Compute the dt
	double previousTime = _solverHandler->getPreviousTime();
	double dt = time - previousTime;

	// Don't do anything if it is not on the stride
	if (((PetscInt)((time + dt / 10.0) / _hdf5Stride) <= _hdf5Previous) &&
		timestep > 0)
		PetscFunctionReturn(0);

	// Update the previous time
	if ((PetscInt)((time + dt / 10.0) / _hdf5Stride) > _hdf5Previous)
		_hdf5Previous++;

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the network and dof
	auto& network = _solverHandler->getNetwork();
	const auto dof = network.getDOF();

	// Create an array for the concentration
	double concArray[dof + 1][2];

	// Open the existing HDF5 file
	auto xolotlComm = util::getMPIComm();
	io::XFile checkpointFile(
		_hdf5OutputName, xolotlComm, io::XFile::AccessMode::OpenReadWrite);

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	CHKERRQ(ierr);

	// Add a concentration time step group for the current time step.
	auto concGroup = checkpointFile.getGroup<io::XFile::ConcentrationGroup>();
	assert(concGroup);
	auto tsGroup = concGroup->addTimestepGroup(
		_loopNumber, timestep, time, previousTime, currentTimeStep);

	// Determine the concentration values we will write.
	io::XFile::TimestepGroup::Concs1DType concs(1);

	// Access the solution data for the current grid point.
	gridPointSolution = solutionArray[0];

	for (auto l = 0; l < dof + 1; ++l) {
		if (std::fabs(gridPointSolution[l]) > 1.0e-16) {
			concs[0].emplace_back(l, gridPointSolution[l]);
		}
	}

	// Write our concentration data to the current timestep group
	// in the HDF5 file.
	// We only write the data for the grid points we own.
	tsGroup->writeConcentrations(checkpointFile, 0, concs);

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::computeXenonRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	using NetworkType = core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Degrees of freedom is the total number of clusters in the network
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	const auto dof = network.getDOF();

	// Get the array of concentration
	PetscReal** solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration and other values over the grid
	double xeConcentration = 0.0, bubbleConcentration = 0.0, radii = 0.0,
		   partialBubbleConcentration = 0.0, partialRadii = 0.0,
		   partialSize = 0.0;

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal* gridPointSolution;

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	// Get the minimum size for the radius
	auto minSizes = _solverHandler->getMinSizes();

	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(gridPointSolution, dof);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof);
	deep_copy(dConcs, hConcs);

	// Get the concentrations
	using TQ = core::network::IReactionNetwork::TotalQuantity;
	using Q = TQ::Type;
	using TQA = util::Array<TQ, 6>;
	auto id = core::network::SpeciesId(Spec::Xe, network.getSpeciesListSize());
	auto ms = static_cast<AmountType>(minSizes[id()]);
	auto totals = network.getTotals(dConcs,
		TQA{TQ{Q::total, id, 1}, TQ{Q::atom, id, 1}, TQ{Q::radius, id, 1},
			TQ{Q::total, id, ms}, TQ{Q::atom, id, ms}, TQ{Q::radius, id, ms}});
	bubbleConcentration = totals[0];
	xeConcentration = totals[1];
	radii = totals[2];
	partialBubbleConcentration = totals[3];
	partialSize = totals[4];
	partialRadii = totals[5];

	// Print the result
	XOLOTL_LOG << "\nTime: " << time << '\n'
			   << "Xenon concentration = " << xeConcentration << "\n\n";

	// Make sure the average partial radius makes sense
	double averagePartialRadius = 0.0, averagePartialSize = 0.0;
	if (partialBubbleConcentration > 1.e-16) {
		averagePartialRadius = partialRadii / partialBubbleConcentration;
		averagePartialSize = partialSize / partialBubbleConcentration;
	}

	// Uncomment to write the content in a file
	std::ofstream outputFile;
	outputFile.open("retentionOut.txt", std::ios::app);
	outputFile << time << " " << xeConcentration << " "
			   << radii / bubbleConcentration << " " << averagePartialRadius
			   << " " << partialBubbleConcentration << " " << averagePartialSize
			   << std::endl;
	outputFile.close();

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::computeHeliumRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	using NetworkType = core::network::IPSIReactionNetwork;

	// Degrees of freedom is the total number of clusters in the network
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	const auto dof = network.getDOF();
	auto numSpecies = network.getSpeciesListSize();

	// Get the array of concentration
	PetscReal** solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration and other values over the grid
	double heConcentration = 0.0, cb = 0.0, avHe = 0.0, avV = 0.0;

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal* gridPointSolution;

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	// Get the minimum size for the radius
	auto minSizes = _solverHandler->getMinSizes();

	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(gridPointSolution, dof);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof);
	deep_copy(dConcs, hConcs);

	// Get the concentrations
	auto id = core::network::SpeciesId(numSpecies);
	heConcentration = network.getTotalAtomConcentration(dConcs, id, 1);
	double vConcentration = network.getTotalAtomConcentration(dConcs, ++id, 1);
	double iConcentration = network.getTotalAtomConcentration(dConcs, ++id, 1);
	cb = gridPointSolution[dof - 3];
	avHe = gridPointSolution[dof - 2];
	avV = gridPointSolution[dof - 1];

	// Print the result
	XOLOTL_LOG << "\nTime: " << time << '\n'
			   << "Helium concentration = " << heConcentration << "\n"
			   << avV / cb << "\n\n";

	// Uncomment to write the content in a file
	std::ofstream outputFile;
	outputFile.open("retentionOut.txt", std::ios::app);
	outputFile << time << " " << heConcentration << " " << vConcentration << " "
			   << iConcentration << " " << cb << " " << avHe / cb << " "
			   << avV / cb << " " << avHe + heConcentration << std::endl;
	outputFile.close();

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::computeAlloy(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the array of concentration
	PetscReal** solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	using NetworkType = core::network::AlloyReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Degrees of freedom is the total number of clusters in the network
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	const auto dof = network.getDOF();
	auto numSpecies = network.getSpeciesListSize();
	auto myData = std::vector<double>(numSpecies * 4, 0.0);

	// Get the minimum size for the loop densities and diameters
	auto minSizes = _solverHandler->getMinSizes();

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal* gridPointSolution;

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(gridPointSolution, dof);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof);
	deep_copy(dConcs, hConcs);

	// Loop on the species
	for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
		using TQ = core::network::IReactionNetwork::TotalQuantity;
		using Q = TQ::Type;
		using TQA = util::Array<TQ, 4>;
		auto ms = static_cast<AmountType>(minSizes[id()]);
		auto totals = network.getTotals(dConcs,
			TQA{TQ{Q::total, id, 1}, TQ{Q::radius, id, 1}, TQ{Q::total, id, ms},
				TQ{Q::radius, id, ms}});

		myData[4 * id()] = totals[0];
		myData[(4 * id()) + 1] = 2.0 * totals[1] / myData[4 * id()];
		myData[(4 * id()) + 2] = totals[2];
		myData[(4 * id()) + 3] = 2.0 * totals[3] / myData[(4 * id()) + 2];
	}

	// Set the output precision
	const int outputPrecision = 5;

	// Open the output file
	std::fstream outputFile;
	outputFile.open("Alloy.dat", std::fstream::out | std::fstream::app);
	outputFile << std::setprecision(outputPrecision);

	// Output the data
	outputFile << timestep << " " << time << " ";
	for (auto i = 0; i < numSpecies; ++i) {
		outputFile << myData[i * 4] << " " << myData[(i * 4) + 1] << " "
				   << myData[(i * 4) + 2] << " " << myData[(i * 4) + 3] << " ";
	}
	outputFile << std::endl;

	// Close the output file
	outputFile.close();

	// Restore the PETSc solution array
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::computeAlphaZr(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the array of concentration
	PetscReal** solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	using NetworkType = core::network::ZrReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Degrees of freedom is the total number of clusters in the network
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	const auto dof = network.getDOF();
	auto numSpecies = network.getSpeciesListSize();
	auto myData = std::vector<double>(numSpecies * 6, 0.0);

	// Get the minimum size for the loop densities and diameters
	auto minSizes = _solverHandler->getMinSizes();

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal* gridPointSolution;

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(gridPointSolution, dof);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof);
	deep_copy(dConcs, hConcs);

	// Loop on the species
	for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
		using TQ = core::network::IReactionNetwork::TotalQuantity;
		using Q = TQ::Type;
		using TQA = util::Array<TQ, 6>;
		auto ms = static_cast<AmountType>(minSizes[id()]);
		auto totals = network.getTotals(dConcs,
			TQA{TQ{Q::total, id, 1}, TQ{Q::atom, id, 1}, TQ{Q::radius, id, 1},
				TQ{Q::total, id, ms}, TQ{Q::atom, id, ms},
				TQ{Q::radius, id, ms}});

		myData[6 * id()] = totals[0];
		myData[6 * id() + 1] = totals[1];
		myData[(6 * id()) + 2] = 2.0 * totals[2] / myData[6 * id()];
		myData[(6 * id()) + 3] = totals[3];
		myData[(6 * id()) + 4] = totals[4];
		myData[(6 * id()) + 5] = 2.0 * totals[5] / myData[(6 * id()) + 3];
	}

	// Set the output precision
	const int outputPrecision = 5;

	// Open the output file
	std::fstream outputFile;
	outputFile.open("AlphaZr.dat", std::fstream::out | std::fstream::app);
	outputFile << std::setprecision(outputPrecision);

	// Output the data
	outputFile << timestep << " " << time << " ";

	for (auto i = 0; i < numSpecies; ++i) {
		outputFile << myData[i * 6] << " " << myData[(i * 6) + 1] << " "
				   << myData[(i * 6) + 2] << " " << myData[(i * 6) + 3] << " "
				   << myData[(i * 6) + 4] << " " << myData[(i * 6) + 5] << " ";
	}

	outputFile << std::endl;

	// Close the output file
	outputFile.close();

	// Restore the PETSc solution array
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::monitorScatter(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % 10 != 0)
		PetscFunctionReturn(0);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the network and its size
	using NetworkType = core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Region = typename NetworkType::Region;
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	auto networkSize = network.getNumClusters();

	// Create a DataPoint vector to store the data to give to the data provider
	// for the visualization
	auto myPoints =
		std::make_shared<std::vector<viz::dataprovider::DataPoint>>();

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	for (auto i = 0; i < networkSize; i++) {
		// Create a DataPoint with the concentration[i] as the value
		// and add it to myPoints
		auto cluster = network.getCluster(i, plsm::HostMemSpace{});
		const Region& clReg = cluster.getRegion();
		for (auto j : makeIntervalRange(clReg[Spec::Xe])) {
			viz::dataprovider::DataPoint aPoint;
			aPoint.value = gridPointSolution[i];
			aPoint.t = time;
			aPoint.x = (double)j;
			myPoints->push_back(aPoint);
		}
	}

	// Get the data provider and give it the points
	_scatterPlot->getDataProvider()->setDataPoints(myPoints);

	// Change the title of the plot and the name of the data
	std::stringstream title;
	title << "Size Distribution";
	_scatterPlot->getDataProvider()->setDataName(title.str());
	_scatterPlot->plotLabelProvider->titleLabel = title.str();
	// Give the time to the label provider
	std::stringstream timeLabel;
	timeLabel << "time: " << std::setprecision(4) << time << "s";
	_scatterPlot->plotLabelProvider->timeLabel = timeLabel.str();
	// Get the current time step
	PetscReal currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	CHKERRQ(ierr);
	// Give the timestep to the label provider
	std::stringstream timeStepLabel;
	timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep << "s";
	_scatterPlot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

	// Render and save in file
	std::stringstream fileName;
	fileName << "Scatter_TS" << timestep << ".png";
	_scatterPlot->render(fileName.str());

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::monitorBubble(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	//	if (timestep % 10 != 0)
	//		PetscFunctionReturn(0);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the network
	using NetworkType = core::network::FeReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	using Region = typename NetworkType::Region;

	// Get the network and its size
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	const auto networkSize = network.getNumClusters();

	// Create the output file
	std::ofstream outputFile;
	std::stringstream name;
	name << "bubble_" << timestep << ".dat";
	outputFile.open(name.str());
	outputFile << "#lo_He hi_He lo_V hi_V conc" << std::endl;

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	// Initialize the total helium and concentration before looping
	double concTot = 0.0, heliumTot = 0.0;

	// Consider each cluster.
	for (auto i = 0; i < networkSize; i++) {
		auto cluster = network.getCluster(i, plsm::HostMemSpace{});
		const Region& clReg = cluster.getRegion();
		Composition lo = clReg.getOrigin();
		Composition hi = clReg.getUpperLimitPoint();

		if (lo.isOnAxis(Spec::I) || lo.isOnAxis(Spec::V) ||
			lo.isOnAxis(Spec::He))
			continue;

		// For compatibility with previous versions, we output
		// the value of a closed upper bound of the He and V intervals.
		outputFile << lo[Spec::He] << " " << hi[Spec::He] - 1 << " "
				   << lo[Spec::V] << " " << hi[Spec::V] - 1 << " "
				   << gridPointSolution[i] << std::endl;
	}

	// Close the file
	outputFile.close();

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
} // namespace monitor
} // namespace solver
} // namespace xolotl
