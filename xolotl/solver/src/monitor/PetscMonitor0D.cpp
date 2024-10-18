#include <petscconf.h>
#include <petscdmda_kokkos.hpp>

#include <xolotl/core/Types.h>
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
#include <xolotl/util/Tokenizer.h>
#include <xolotl/viz/dataprovider/CvsXDataProvider.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
PetscErrorCode
monitorBubble(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<PetscMonitor0D*>(ictx)->monitorBubble(
		ts, timestep, time, solution));
	PetscFunctionReturn(0);
}

void
PetscMonitor0D::setup(int loop)
{
	_loopNumber = loop;

	// Get xolotlViz handler registry
	auto vizHandlerRegistry = _solverHandler->getVizHandler();

	// Flags to launch the monitors or not
	PetscBool flagCheck, flag1DPlot, flagBubble, flagStatus, flagAlloy,
		flagXeRetention, flagHeRetention, flagLargest, flagZr;

	// Check the option -check_collapse
	PetscCallVoid(
		PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCheck));

	// Check the option -plot_1d
	PetscCallVoid(PetscOptionsHasName(NULL, NULL, "-plot_1d", &flag1DPlot));

	// Check the option -start_stop
	PetscCallVoid(PetscOptionsHasName(NULL, NULL, "-start_stop", &flagStatus));

	// Check the option -bubble
	PetscCallVoid(PetscOptionsHasName(NULL, NULL, "-bubble", &flagBubble));

	// Check the option -alloy
	PetscCallVoid(PetscOptionsHasName(NULL, NULL, "-alloy", &flagAlloy));

	// Check the option -alpha_zr
	PetscCallVoid(PetscOptionsHasName(NULL, NULL, "-alpha_zr", &flagZr));

	// Check the option -xenon_retention
	PetscCallVoid(
		PetscOptionsHasName(NULL, NULL, "-xenon_retention", &flagXeRetention));

	// Check the option -helium_retention
	PetscCallVoid(
		PetscOptionsHasName(NULL, NULL, "-helium_retention", &flagHeRetention));

	// Check the option -largest_conc
	PetscCallVoid(
		PetscOptionsHasName(NULL, NULL, "-largest_conc", &flagLargest));

	// Determine if we have an existing restart file,
	// and if so, it it has had timesteps written to it.
	std::unique_ptr<io::XFile> networkFile;
	std::unique_ptr<io::XFile::TimestepGroup> lastTsGroup;
	bool hasConcentrations = false;
	if (_solverHandler->checkForRestart()) {
		auto restartFilePath = _solverHandler->getRestartFilePath();
		networkFile = std::make_unique<io::XFile>(restartFilePath);
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
		PetscCallVoid(PetscOptionsGetReal(
			NULL, NULL, "-check_collapse", &timeStepThreshold, &flag));
		if (!flag)
			timeStepThreshold = 1.0e-16;

		// Set the post step process that tells the solver when to stop if the
		// time step collapse
		PetscCallVoid(TSSetPostStep(_ts, checkTimeStep));
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
		PetscCallVoid(
			TSMonitorSet(_ts, monitor::monitorScatter, this, nullptr));
	}

	// Set the monitor to save text file of the mean concentration of bubbles
	if (flagBubble) {
		// monitorBubble0D will be called at each timestep
		PetscCallVoid(TSMonitorSet(_ts, monitor::monitorBubble, this, nullptr));
	}

	// Set the monitor to output data for Alloy
	if (flagAlloy) {
		_solverHandler->getNetwork().writeMonitorOutputHeader();

		// computeAlloy0D will be called at each timestep
		PetscCallVoid(TSMonitorSet(_ts, monitor::computeAlloy, this, nullptr));
	}
	// Set the monitor to output data for AlphaZr
	if (flagZr) {
		_solverHandler->getNetwork().writeMonitorOutputHeader();

		// computeAlphaZr will be called at each timestep
		PetscCallVoid(
			TSMonitorSet(_ts, monitor::computeAlphaZr, this, nullptr));
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
			auto fluences = lastTsGroup->readFluence();
			fluxHandler->setFluence(fluences);
		}

		// computeXenonRetention0D will be called at each timestep
		PetscCallVoid(
			TSMonitorSet(_ts, monitor::computeXenonRetention, this, nullptr));

		using NetworkType = core::network::NEReactionNetwork;
		using Spec = typename NetworkType::Species;
		using Composition = typename NetworkType::Composition;
		using Region = typename NetworkType::Region;
		auto& network =
			dynamic_cast<NetworkType&>(_solverHandler->getNetwork());

		// Uncomment to clear the file where the retention will be written
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt");
		outputFile << "#time content ";

		std::ifstream reactionFile;
		reactionFile.open(_solverHandler->getReactionFilePath());
		// Get the line
		std::string line;
		getline(reactionFile, line);
		// Read the first line
		std::vector<double> tokens;
		util::Tokenizer<double>{line}(tokens);
		// And start looping on the lines
		while (tokens.size() > 0) {
			// Find the Id of the cluster
			IdType nXe = static_cast<IdType>(tokens[0]);
			IdType nV = static_cast<IdType>(tokens[1]);
			IdType nI = static_cast<IdType>(tokens[2]);
			auto comp =
				std::vector<AmountType>(network.getSpeciesListSize(), 0);
			auto clusterSpecies = network.parseSpeciesId("Xe");
			comp[clusterSpecies()] = nXe;
			clusterSpecies = network.parseSpeciesId("V");
			comp[clusterSpecies()] = nV;
			clusterSpecies = network.parseSpeciesId("I");
			comp[clusterSpecies()] = nI;

			auto clusterId = network.findClusterId(comp);
			// Check that it is present in the network
			if (clusterId != NetworkType::invalidIndex()) {
				_clusterOrder.push_back(clusterId);
				if (nI > 0)
					outputFile << "I_" << nI << " ";
				else if (nV > 0 and nXe == 0)
					outputFile << "V_" << nV << " ";
				else if (nXe > 0 and nV == 0)
					outputFile << "Xe_" << nXe << " ";
				else
					outputFile << "Xe_" << nXe << "V_" << nV << " ";
			}

			getline(reactionFile, line);
			if (line == "Reactions")
				break;

			tokens = util::Tokenizer<double>{line}();
		}
		outputFile << "Xe/SD var" << std::endl;
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
		PetscCallVoid(
			TSMonitorSet(_ts, monitor::computeHeliumRetention, this, nullptr));

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
		PetscCallVoid(PetscOptionsGetReal(
			NULL, NULL, "-largest_conc", &_largestThreshold, &flag));

		// monitorLargest1D will be called at each timestep
		PetscCallVoid(
			TSMonitorSet(_ts, monitor::monitorLargest, this, nullptr));
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		PetscCallVoid(PetscOptionsGetReal(
			NULL, NULL, "-start_stop", &_hdf5Stride, &flag));
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
		if (_hdf5OutputName != _solverHandler->getRestartFilePath()) {
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
			writeNetwork(xolotlComm, _hdf5OutputName,
				_solverHandler->getRestartFilePath());
		}

		// startStop0D will be called at each timestep
		PetscCallVoid(TSMonitorSet(_ts, monitor::startStop, this, nullptr));
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	PetscCallVoid(TSMonitorSet(_ts, monitor::monitorTime, this, nullptr));
}

PetscErrorCode
PetscMonitor0D::monitorLargest(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declaration
	double **solutionArray, *gridPointSolution;

	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	PetscCall(TSGetDM(ts, &da));

	// Get the solutionArray
	PetscCall(DMDAVecGetArrayDOF(da, solution, &solutionArray));

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];
	// Check the concentration
	if (gridPointSolution[_largestClusterId] > _largestThreshold) {
		PetscCall(TSSetConvergedReason(ts, TS_CONVERGED_USER));
		// Send an error
		throw std::runtime_error(
			"\nxolotlSolver::Monitor0D: The largest cluster "
			"concentration is too high!!");
	}

	// Restore the solutionArray
	PetscCall(DMDAVecRestoreArrayDOF(da, solution, &solutionArray));

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::startStopImpl(TS ts, PetscInt timestep, PetscReal time,
	Vec solution, io::XFile& checkpointFile, io::XFile::TimestepGroup* tsGroup,
	[[maybe_unused]] const std::vector<std::string>& speciesNames)
{
	// Initial declaration
	const double **solutionArray, *gridPointSolution;

	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	PetscCall(TSGetDM(ts, &da));

	// Get the solutionArray
	PetscCall(DMDAVecGetArrayDOFRead(da, solution, &solutionArray));

	// Get the network and dof
	auto& network = _solverHandler->getNetwork();
	const auto dof = network.getDOF();

	// Create an array for the concentration
	double concArray[dof + 1][2];

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
	PetscCall(DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray));

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::computeXenonRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	PetscCall(TSGetDM(ts, &da));

	using NetworkType = core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Degrees of freedom is the total number of clusters in the network
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());

	// Get the array of concentration
	PetscOffsetView<const PetscReal**> solutionArray;
	PetscCall(DMDAVecGetKokkosOffsetViewDOF(da, solution, &solutionArray));

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal* gridPointSolution;
	PetscReal** solutionArrayH;
	PetscCall(DMDAVecGetArrayDOFRead(da, solution, &solutionArrayH));
	gridPointSolution = solutionArrayH[0];

	// Store the concentration and other values over the grid
	double xeConcentration = 0.0;

	// Get the pointer to the beginning of the solution data for this grid point
	auto concs = subview(solutionArray, 0, Kokkos::ALL).view();

	// Get the minimum size for the radius
	auto minSizes = _solverHandler->getMinSizes();

	// Get the concentrations
	xeConcentration = network.getTotalAtomConcentration(concs, Spec::Xe, 1);

	// Print the result
	XOLOTL_LOG << "\nTime: " << time << '\n'
			   << "Xenon concentration = " << xeConcentration << "\n\n";

	// Uncomment to write the content in a file
	constexpr double k_B = ::xolotl::core::kBoltzmann;
	std::ofstream outputFile;
	outputFile.open("retentionOut.txt", std::ios::app);
	outputFile << time << " " << xeConcentration << " ";
	for (auto id : _clusterOrder) {
		outputFile << gridPointSolution[id] << " ";
	}
	if (xeConcentration < 1.0e-16)
		outputFile << "0 0" << std::endl;
	else {
		auto ratio = network.getTotalVolumeRatio(concs, Spec::Xe, 2);
		auto variance =
			network.getTotalRatioVariance(concs, Spec::Xe, ratio, 2);
		outputFile << ratio << " " << variance << std::endl;
	}
	outputFile.close();

	// Restore the solutionArray
	PetscCall(DMDAVecRestoreKokkosOffsetViewDOF(da, solution, &solutionArray));

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
	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	PetscCall(TSGetDM(ts, &da));

	// Get the array of concentration
	PetscOffsetView<const PetscReal**> concs;
	PetscCall(DMDAVecGetKokkosOffsetViewDOF(da, solution, &concs));
	auto concOffset = subview(concs, 0, Kokkos::ALL).view();

	using NetworkType = core::network::AlloyReactionNetwork;
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());

	auto myData = network.getMonitorDataValues(concOffset, 1.0);

	network.writeMonitorDataLine(myData, time);

	// Restore the PETSc solution array
	PetscCall(DMDAVecRestoreKokkosOffsetViewDOF(da, solution, &concs));

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::computeAlphaZr(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	PetscCall(TSGetDM(ts, &da));

	// Get the array of concentration
	PetscOffsetView<const PetscScalar**> concs;
	PetscCall(DMDAVecGetKokkosOffsetViewDOF(da, solution, &concs));
	auto concOffset = subview(concs, 0, Kokkos::ALL).view();

	using NetworkType = core::network::ZrReactionNetwork;
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());

	auto myData = network.getMonitorDataValues(concOffset, 1.0);

	network.writeMonitorDataLine(myData, time);

	// Restore the PETSc solution array
	PetscCall(DMDAVecRestoreKokkosOffsetViewDOF(da, solution, &concs));

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::monitorScatter(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	double **solutionArray, *gridPointSolution;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % 10 != 0)
		PetscFunctionReturn(0);

	// Get the da from ts
	DM da;
	PetscCall(TSGetDM(ts, &da));

	// Get the solutionArray
	PetscCall(DMDAVecGetArrayDOFRead(da, solution, &solutionArray));

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
	PetscCall(TSGetTimeStep(ts, &currentTimeStep));
	// Give the timestep to the label provider
	std::stringstream timeStepLabel;
	timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep << "s";
	_scatterPlot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

	// Render and save in file
	std::stringstream fileName;
	fileName << "Scatter_TS" << timestep << ".png";
	_scatterPlot->render(fileName.str());

	// Restore the solutionArray
	PetscCall(DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray));

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor0D::monitorBubble(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declaration
	double **solutionArray, *gridPointSolution;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	//	if (timestep % 10 != 0)
	//		PetscFunctionReturn(0);

	// Get the da from ts
	DM da;
	PetscCall(TSGetDM(ts, &da));

	// Get the solutionArray
	PetscCall(DMDAVecGetArrayDOFRead(da, solution, &solutionArray));

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
	PetscCall(DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray));

	PetscFunctionReturn(0);
}
} // namespace monitor
} // namespace solver
} // namespace xolotl
