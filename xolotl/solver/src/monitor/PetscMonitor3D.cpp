#include <xolotl/core/Constants.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/io/XFile.h>
#include <xolotl/perf/ScopedTimer.h>
#include <xolotl/solver/PetscSolver.h>
#include <xolotl/solver/monitor/PetscMonitor3D.h>
#include <xolotl/solver/monitor/PetscMonitorFunctions.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/viz/dataprovider/CvsXDataProvider.h>
#include <xolotl/viz/dataprovider/CvsXYDataProvider.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
PetscErrorCode
monitorSurfaceXY(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<PetscMonitor3D*>(ictx)->monitorSurfaceXY(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
monitorSurfaceXZ(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<PetscMonitor3D*>(ictx)->monitorSurfaceXZ(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

void
PetscMonitor3D::setup()
{
	PetscErrorCode ierr;

	// Get the process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the xolotlViz handler registry
	auto vizHandlerRegistry = _solverHandler->getVizHandler();

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

	// Get the network and its size
	auto& network = _solverHandler->getNetwork();

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

	// Get the da from _ts
	DM da;
	ierr = TSGetDM(_ts, &da);
	checkPetscError(ierr, "setupPetsc3DMonitor: TSGetDM failed.");

	// Get the total size of the grid
	PetscInt Mx, My, Mz;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	checkPetscError(ierr, "setupPetsc3DMonitor: DMDAGetInfo failed.");

	// Get the number of species
	auto numSpecies = network.getSpeciesListSize();

	// Create data depending on the boundary conditions
	if (_solverHandler->getLeftOffset() == 1) {
		_nSurf = std::vector<std::vector<std::vector<double>>>(numSpecies,
			std::vector<std::vector<double>>(My, std::vector<double>(Mz, 0.0)));
		_previousSurfFlux = std::vector<std::vector<std::vector<double>>>(
			numSpecies,
			std::vector<std::vector<double>>(My, std::vector<double>(Mz, 0.0)));
	}
	if (_solverHandler->getRightOffset() == 1) {
		_nBulk = std::vector<std::vector<std::vector<double>>>(numSpecies,
			std::vector<std::vector<double>>(My, std::vector<double>(Mz, 0.0)));
		_previousBulkFlux = std::vector<std::vector<std::vector<double>>>(
			numSpecies,
			std::vector<std::vector<double>>(My, std::vector<double>(Mz, 0.0)));
	}

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
		ierr = TSSetPostStep(_ts, checkTimeStep);
		checkPetscError(
			ierr, "setupPetsc3DMonitor: TSSetPostStep (checkTimeStep) failed.");
	}

	// If the user wants the surface to be able to move
	if (_solverHandler->moveSurface()) {
		// Clear the vector just in case
		_iClusterIds.clear();

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
				_iClusterIds.push_back(clusterId);
				iSize++;
			}
			else
				iClusterExists = false;
		}

		// Get the interstitial information at the surface if concentrations
		// were stored
		if (hasConcentrations) {
			assert(lastTsGroup);

			// Get the names of the species in the network
			std::vector<std::string> names;
			for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
				names.push_back(network.getSpeciesName(id));
			}

			// Loop on the names
			for (auto i = 0; i < names.size(); i++) {
				// Create the n attribute name
				std::ostringstream nName;
				nName << "n" << names[i] << "Surf";
				// Read quantity attribute
				_nSurf[i] = lastTsGroup->readData3D(nName.str());

				// Create the previous flux attribute name
				std::ostringstream prevFluxName;
				prevFluxName << "previousFlux" << names[i] << "Surf";
				// Read the attribute
				_previousSurfFlux[i] =
					lastTsGroup->readData3D(prevFluxName.str());
			}

			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			_solverHandler->setPreviousTime(previousTime);
		}

		// Get the sputtering yield
		_sputteringYield = _solverHandler->getSputteringYield();

		// Master process
		if (procId == 0) {
			// Clear the file where the surface will be written
			std::ofstream outputFile;
			outputFile.open("surface.txt");
			outputFile << "#time heights" << std::endl;
			outputFile.close();
		}

		// Set directions and terminate flags for the surface event
		PetscInt direction[1];
		PetscBool terminate[1];
		direction[0] = 0;
		terminate[0] = PETSC_FALSE;
		// Set the TSEvent
		ierr = TSSetEventHandler(_ts, 1, direction, terminate,
			monitor::eventFunction, monitor::postEventFunction, this);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: TSSetEventHandler (eventFunction3D) failed.");
	}

	// Set the monitor to save performance plots (has to be in parallel)
	if (flagPerf) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a ScatterPlot
			_perfPlot = vizHandlerRegistry->getPlot(viz::PlotType::SCATTER);

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>();
			labelProvider->axis1Label = "Process ID";
			labelProvider->axis2Label = "Solver Time";

			// Give it to the plot
			_perfPlot->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider =
				std::make_shared<viz::dataprovider::CvsXDataProvider>();

			// Give it to the plot
			_perfPlot->setDataProvider(dataProvider);
		}

		// monitorPerf will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorPerf, this, nullptr);
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
			_solverHandler->setPreviousTime(previousTime);
			// Initialize the fluence
			auto fluxHandler = _solverHandler->getFluxHandler();
			// Increment the fluence with the value at this current timestep
			fluxHandler->computeFluence(previousTime);

			// Get the names of the species in the network
			std::vector<std::string> names;
			for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
				names.push_back(network.getSpeciesName(id));
			}

			// If the surface is a free surface
			if (_solverHandler->getLeftOffset() == 1) {
				// Loop on the names
				for (auto i = 0; i < names.size(); i++) {
					// Create the n attribute name
					std::ostringstream nName;
					nName << "n" << names[i] << "Surf";
					// Read quantity attribute
					_nSurf[i] = lastTsGroup->readData3D(nName.str());

					// Create the previous flux attribute name
					std::ostringstream prevFluxName;
					prevFluxName << "previousFlux" << names[i] << "Surf";
					// Read the attribute
					_previousSurfFlux[i] =
						lastTsGroup->readData3D(prevFluxName.str());
				}
			}

			// If the bottom is a free surface
			if (_solverHandler->getRightOffset() == 1) {
				// Loop on the names
				for (auto i = 0; i < names.size(); i++) {
					// Create the n attribute name
					std::ostringstream nName;
					nName << "n" << names[i] << "Bulk";
					// Read quantity attribute
					_nBulk[i] = lastTsGroup->readData3D(nName.str());

					// Create the previous flux attribute name
					std::ostringstream prevFluxName;
					prevFluxName << "previousFlux" << names[i] << "Bulk";
					// Read the attribute
					_previousBulkFlux[i] =
						lastTsGroup->readData3D(prevFluxName.str());
				}
			}
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeFluence, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc3DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeHeliumRetention3D will be called at each timestep
		ierr =
			TSMonitorSet(_ts, monitor::computeHeliumRetention, this, nullptr);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: TSMonitorSet (computeHeliumRetention) "
			"failed.");

		// Master process
		if (procId == 0) {
			auto numSpecies = network.getSpeciesListSize();
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("retentionOut.txt");
			outputFile << "#time fluence ";
			for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
				auto speciesName = network.getSpeciesName(id);
				outputFile << speciesName << "_content ";
			}
			if (_solverHandler->getRightOffset() == 1) {
				for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
					auto speciesName = network.getSpeciesName(id);
					outputFile << speciesName << "_bulk ";
				}
			}
			if (_solverHandler->getLeftOffset() == 1) {
				for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
					auto speciesName = network.getSpeciesName(id);
					outputFile << speciesName << "_surface ";
				}
			}
			outputFile << std::endl;
			outputFile.close();
		}
	}

	// Set the monitor to compute the xenon fluence and the retention
	// for the retention calculation
	if (flagXeRetention) {
		// Get the da from _ts
		DM da;
		ierr = TSGetDM(_ts, &da);
		checkPetscError(ierr, "setupPetsc3DMonitor: TSGetDM failed.");
		// Get the local boundaries
		PetscInt xm, ym, zm;
		ierr = DMDAGetCorners(da, NULL, NULL, NULL, &xm, &ym, &zm);
		checkPetscError(ierr, "setupPetsc3DMonitor: DMDAGetCorners failed.");
		// Create the local vectors on each process
		_solverHandler->createLocalNE(xm, ym, zm);

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

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeFluence, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc3DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeXenonRetention3D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeXenonRetention, this, nullptr);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: TSMonitorSet (computeXenonRetention) "
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
			_surfacePlotXY =
				vizHandlerRegistry->getPlot(viz::PlotType::SURFACE);

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>();
			labelProvider->axis1Label = "Depth (nm)";
			labelProvider->axis2Label = "Y (nm)";
			labelProvider->axis3Label = "Concentration";

			// Give it to the plot
			_surfacePlotXY->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider =
				std::make_shared<viz::dataprovider::CvsXYDataProvider>();

			// Give it to the plot
			_surfacePlotXY->setDataProvider(dataProvider);
		}

		// monitorSurfaceXY will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorSurfaceXY, this, nullptr);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: TSMonitorSet (monitorSurfaceXY) failed.");
	}

	// Set the monitor to save surface plots of clusters concentration
	if (flag2DXZPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a SurfacePlot
			_surfacePlotXZ =
				vizHandlerRegistry->getPlot(viz::PlotType::SURFACE);

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>();
			labelProvider->axis1Label = "Depth (nm)";
			labelProvider->axis2Label = "Z (nm)";
			labelProvider->axis3Label = "Concentration";

			// Give it to the plot
			_surfacePlotXZ->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider =
				std::make_shared<viz::dataprovider::CvsXYDataProvider>();

			// Give it to the plot
			_surfacePlotXZ->setDataProvider(dataProvider);
		}

		// monitorSurfaceXZ will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorSurfaceXZ, this, nullptr);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: TSMonitorSet (monitorSurfaceXZ) failed.");
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
			"setupPetsc3DMonitor: PetscOptionsGetReal (-largest_conc) failed.");

		// monitorLargest3D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorLargest, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc3DMonitor: TSMonitorSet (monitorLargest) failed.");
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr =
			PetscOptionsGetReal(NULL, NULL, "-start_stop", &_hdf5Stride, &flag);
		checkPetscError(ierr,
			"setupPetsc3DMonitor: PetscOptionsGetReal (-start_stop) failed.");
		if (!flag)
			_hdf5Stride = 1.0;

		// Compute the correct _hdf5Previous for a restart
		if (hasConcentrations) {
			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			_solverHandler->setPreviousTime(previousTime);
			_hdf5Previous = (PetscInt)(previousTime / _hdf5Stride);
		}

		// Don't do anything if both files have the same name
		if (_hdf5OutputName != _solverHandler->getNetworkName()) {
			// Get the physical grid in the x direction
			auto grid = _solverHandler->getXGrid();

			// Setup step size variables
			double hy = _solverHandler->getStepSizeY();
			double hz = _solverHandler->getStepSizeZ();

			// Create a checkpoint file.
			// Create and initialize a checkpoint file.
			// We do this in its own scope so that the file
			// is closed when the file object goes out of scope.
			// We want it to close before we (potentially) copy
			// the network from another file using a single-process
			// MPI communicator.
			{
				io::XFile checkpointFile(
					_hdf5OutputName, grid, xolotlComm, My, hy, Mz, hz);
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

		// startStop3D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::startStop, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc3DMonitor: TSMonitorSet (startStop) failed.");
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(_ts, monitor::monitorTime, this, nullptr);
	checkPetscError(
		ierr, "setupPetsc3DMonitor: TSMonitorSet (monitorTime) failed.");
}

PetscErrorCode
PetscMonitor3D::monitorLargest(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
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

	// Get the local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Loop on the local grid
	for (auto k = zs; k < zs + zm; k++)
		for (auto j = ys; j < ys + ym; j++)
			for (auto i = xs; i < xs + xm; i++) {
				// Get the pointer to the beginning of the solution data for
				// this grid point
				gridPointSolution = solutionArray[k][j][i];
				// Check the concentration
				if (gridPointSolution[_largestClusterId] > _largestThreshold) {
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

PetscErrorCode
PetscMonitor3D::startStop(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	const double ****solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Get the local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

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
	auto& network = _solverHandler->getNetwork();
	const auto dof = network.getDOF();

	// Create an array for the concentration
	double concArray[dof + 1][2];

	// Get the vector of positions of the surface
	std::vector<std::vector<int>> surfaceIndices;
	for (auto i = 0; i < My; i++) {
		// Create a temporary vector
		std::vector<int> temp;
		for (auto j = 0; j < Mz; j++) {
			temp.push_back(_solverHandler->getSurfacePosition(i, j));
		}
		// Add the temporary vector to the vector of surface indices
		surfaceIndices.push_back(temp);
	}

	// Open the existing HDF5 file.
	io::XFile checkpointFile(
		_hdf5OutputName, xolotlComm, io::XFile::AccessMode::OpenReadWrite);

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	CHKERRQ(ierr);

	// Add a concentration sub group
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

	if (_solverHandler->moveSurface() || _solverHandler->getLeftOffset() == 1) {
		// Write the surface positions and the associated interstitial
		// quantities in the concentration sub group
		tsGroup->writeSurface3D(
			surfaceIndices, _nSurf, _previousSurfFlux, names);
	}

	// Write the bottom impurity information if the bottom is a free surface
	if (_solverHandler->getRightOffset() == 1)
		tsGroup->writeBottom3D(_nBulk, _previousBulkFlux, names);

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

PetscErrorCode
PetscMonitor3D::computeHeliumRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Get the local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = _solverHandler->getFluxHandler();
	// Get the diffusion handler
	auto diffusionHandler = _solverHandler->getDiffusionHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the physical grid in the x direction
	auto grid = _solverHandler->getXGrid();

	// Get the network
	using NetworkType = core::network::IPSIReactionNetwork;
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	const auto dof = network.getDOF();

	// Setup step size variables
	double hy = _solverHandler->getStepSizeY();
	double hz = _solverHandler->getStepSizeZ();

	// Get the array of concentration
	double ****solutionArray, *gridPointSolution;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration over the grid
	auto numSpecies = network.getSpeciesListSize();
	auto specIdI = network.getInterstitialSpeciesId();
	auto myConcData = std::vector<double>(numSpecies, 0.0);

	// Loop on the grid
	for (auto zk = zs; zk < zs + zm; zk++) {
		for (auto yj = ys; yj < ys + ym; yj++) {
			// Get the surface position
			auto surfacePos = _solverHandler->getSurfacePosition(yj, zk);

			for (auto xi = xs; xi < xs + xm; xi++) {
				// Boundary conditions
				if (xi < surfacePos + _solverHandler->getLeftOffset() ||
					xi >= Mx - _solverHandler->getRightOffset())
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

	// Get the delta time from the previous timestep to this timestep
	double previousTime = _solverHandler->getPreviousTime();
	double dt = time - previousTime;

	// Look at the fluxes going in the bulk if the bottom is a free surface
	if (_solverHandler->getLeftOffset() == 1) {
		// Get the vector of diffusing clusters
		auto diffusingIds = diffusionHandler->getDiffusingIds();

		// Loop on every Y, Z position
		for (auto k = 0; k < Mz; k++)
			for (auto j = 0; j < My; j++) {
				// Set the surface position
				auto xi = _solverHandler->getSurfacePosition(j, k) + 1;
				// Value to know on which processor is the bottom
				int surfProc = 0;

				// Check we are on the right proc
				if (xi >= xs && xi < xs + xm && j >= ys && j < ys + ym &&
					k >= zs && k < zs + zm) {
					// Compute the total number of impurities that left at the
					// surface
					if (timestep > 0) {
						for (auto i = 0; i < numSpecies; ++i) {
							if (i == specIdI() && _solverHandler->moveSurface())
								continue;
							_nSurf[i][j][k] += _previousSurfFlux[i][j][k] * dt;
						}
					}
					auto myFluxData = std::vector<double>(numSpecies, 0.0);

					// Get the pointer to the beginning of the solution data for
					// this grid point
					gridPointSolution = solutionArray[k][j][xi];

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
					double factor = 2.0 * hy * hz / (hxLeft + hxRight);

					network.updateOutgoingDiffFluxes(gridPointSolution, factor,
						diffusingIds, myFluxData, xi - xs);

					// Take into account the surface advection
					// Get the surface advection handler
					auto advecHandler = _solverHandler->getAdvectionHandler();
					// Get the sink strengths and advecting clusters
					auto sinkStrengths = advecHandler->getSinkStrengths();
					auto advecClusters = advecHandler->getAdvectingClusters();
					// Set the distance from the surface
					double distance = (grid[xi] + grid[xi + 1]) / 2.0 -
						grid[1] - advecHandler->getLocation();

					network.updateOutgoingAdvecFluxes(gridPointSolution,
						3.0 * hy * hz /
							(core::kBoltzmann * distance * distance * distance *
								distance),
						advecClusters, sinkStrengths, myFluxData, xi - xs);

					for (auto i = 0; i < numSpecies; ++i) {
						if (i == specIdI() && _solverHandler->moveSurface())
							continue;
						_previousSurfFlux[i][j][k] = myFluxData[i];
					}

					// Set the bottom processor
					surfProc = procId;
				}

				// Get which processor will send the information
				int surfId = 0;
				MPI_Allreduce(
					&surfProc, &surfId, 1, MPI_INT, MPI_SUM, xolotlComm);

				// Send the information about impurities
				// to the other processes
				std::vector<double> countFluxData;
				for (auto i = 0; i < numSpecies; ++i) {
					countFluxData.push_back(_nSurf[i][j][k]);
					countFluxData.push_back(_previousSurfFlux[i][j][k]);
				}
				MPI_Bcast(countFluxData.data(), countFluxData.size(),
					MPI_DOUBLE, surfId, xolotlComm);

				// Extract impurity data from broadcast buffer.
				for (auto i = 0; i < numSpecies; ++i) {
					_nSurf[i][j][k] = countFluxData[2 * i];
					_previousSurfFlux[i][j][k] = countFluxData[(2 * i) + 1];
				}
			}
	}

	// Look at the fluxes going in the bulk if the bottom is a free surface
	if (_solverHandler->getRightOffset() == 1) {
		// Set the bottom surface position
		auto xi = Mx - 2;

		// Get the vector of diffusing clusters
		auto diffusingIds = diffusionHandler->getDiffusingIds();

		// Loop on every Z, Y position
		for (auto k = 0; k < Mz; k++)
			for (auto j = 0; j < My; j++) {
				// Value to know on which processor is the bottom
				int bottomProc = 0;

				// Check we are on the right proc
				if (xi >= xs && xi < xs + xm && j >= ys && j < ys + ym &&
					k >= zs && k < zs + zm) {
					// Compute the total number of impurities that left at the
					// surface
					if (timestep > 0) {
						for (auto i = 0; i < numSpecies; ++i) {
							_nBulk[i][j][k] += _previousBulkFlux[i][j][k] * dt;
						}
					}
					auto myFluxData = std::vector<double>(numSpecies, 0.0);

					// Get the pointer to the beginning of the solution data for
					// this grid point
					gridPointSolution = solutionArray[k][j][xi];

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
					double factor = 2.0 * hy * hz / (hxLeft + hxRight);

					network.updateOutgoingDiffFluxes(gridPointSolution, factor,
						diffusingIds, myFluxData, xi - xs);

					for (auto i = 0; i < numSpecies; ++i) {
						_previousBulkFlux[i][j][k] = myFluxData[i];
					}

					// Set the bottom processor
					bottomProc = procId;
				}

				// Get which processor will send the information
				int bottomId = 0;
				MPI_Allreduce(
					&bottomProc, &bottomId, 1, MPI_INT, MPI_SUM, xolotlComm);

				// Send the information about impurities
				// to the other processes
				std::vector<double> countFluxData;
				for (auto i = 0; i < numSpecies; ++i) {
					countFluxData.push_back(_nBulk[i][j][k]);
					countFluxData.push_back(_previousBulkFlux[i][j][k]);
				}
				MPI_Bcast(countFluxData.data(), countFluxData.size(),
					MPI_DOUBLE, bottomId, xolotlComm);

				// Extract impurity data from broadcast buffer.
				for (auto i = 0; i < numSpecies; ++i) {
					_nBulk[i][j][k] = countFluxData[2 * i];
					_previousBulkFlux[i][j][k] = countFluxData[(2 * i) + 1];
				}
			}
	}

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
		auto totalBulk = std::vector<double>(numSpecies, 0.0);
		// Look if the bottom is a free surface
		if (_solverHandler->getRightOffset() == 1) {
			for (auto i = 0; i < numSpecies; ++i) {
				for (auto k = 0; k < Mz; k++)
					for (auto j = 0; j < My; j++) {
						totalBulk[i] += _nBulk[i][j][k];
					}
				totalBulk[i] = totalBulk[i] / surface;
			}
		}
		auto totalSurf = std::vector<double>(numSpecies, 0.0);
		// Look if the bottom is a free surface
		if (_solverHandler->getLeftOffset() == 1) {
			for (auto i = 0; i < numSpecies; ++i) {
				for (auto k = 0; k < Mz; k++)
					for (auto j = 0; j < My; j++) {
						totalSurf[i] += _nSurf[i][j][k];
					}
				totalSurf[i] = totalSurf[i] / surface;
			}
		}

		// Get the fluence
		double fluence = fluxHandler->getFluence();

		// Print the result
		util::StringStream ss;
		ss << "\nTime: " << time << '\n';
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			ss << network.getSpeciesName(id)
			   << " content = " << totalConcData[id()] << '\n';
		}
		ss << "Fluence = " << fluence << "\n\n";
		XOLOTL_LOG << ss.str();

		// Uncomment to write the retention and the fluence in a file
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt", std::ios::app);
		outputFile << time << ' ' << fluence << " ";
		for (auto i = 0; i < numSpecies; ++i) {
			outputFile << totalConcData[i] << " ";
		}
		if (_solverHandler->getRightOffset() == 1) {
			for (auto i = 0; i < numSpecies; ++i) {
				outputFile << totalBulk[i] << ' ';
			}
		}
		if (_solverHandler->getLeftOffset() == 1) {
			for (auto i = 0; i < numSpecies; ++i) {
				outputFile << totalSurf[i] << ' ';
			}
		}
		outputFile << std::endl;
		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor3D::computeXenonRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();

	// Setup step size variables
	double hy = _solverHandler->getStepSizeY();
	double hz = _solverHandler->getStepSizeZ();

	using NetworkType = core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Degrees of freedom is the total number of clusters in the network
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
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
	auto minSizes = _solverHandler->getMinSizes();

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
				_solverHandler->setVolumeFraction(
					volumeFrac, xi - xs, yj - ys, zk - zs);
				_solverHandler->setMonomerConc(
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
	double dt = time - _solverHandler->getPreviousTime();
	// Sum and gather the previous flux
	double globalXeFlux = 0.0;
	// Get the vector from the solver handler
	auto gbVector = _solverHandler->getGBVector();
	// Get the previous flux vector
	auto& localNE = _solverHandler->getLocalNE();
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
			_solverHandler->setLocalXeRate(
				previousXeFlux * dt, xi - xs, yj - ys, zk - zs);
		}
	}
	double totalXeFlux = 0.0;
	MPI_Reduce(
		&globalXeFlux, &totalXeFlux, 1, MPI_DOUBLE, MPI_SUM, 0, xolotlComm);
	// Master process
	if (procId == 0) {
		// Get the previous value of Xe that went to the GB
		double nXenon = _solverHandler->getNXeGB();
		// Compute the total number of Xe that went to the GB
		nXenon += totalXeFlux * dt;
		_solverHandler->setNXeGB(nXenon);
	}

	// Get the number of species
	auto numSpecies = network.getSpeciesListSize();

	// Get the vector of diffusing clusters
	auto diffusionHandler = _solverHandler->getDiffusionHandler();
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
			_solverHandler->setPreviousXeFlux(
				myRate[0], xi - xs, yj - ys, zk - zs);
		}
	}

	// Master process
	if (procId == 0) {
		// Compute the total surface irradiated
		double surface = (double)My * hy * (double)Mz * hz;
		// Get the number of Xe that went to the GB
		double nXenon = _solverHandler->getNXeGB();

		totalConcData[0] = totalConcData[0] / surface;

		// Print the result
		XOLOTL_LOG << "\nTime: " << time << '\n'
				   << "Xenon concentration = " << totalConcData[0] << '\n'
				   << "Xenon GB = " << nXenon / surface << "\n\n";

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

PetscErrorCode
PetscMonitor3D::eventFunction(
	TS ts, PetscReal time, Vec solution, PetscScalar* fvalue)
{
	// Initial declaration
	PetscErrorCode ierr;
	double ****solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	fvalue[0] = 1.0;
	_depthPositions.clear();

	PetscInt TSNumber = -1;
	ierr = TSGetStepNumber(ts, &TSNumber);

	// Skip if it is the same TS as before
	if (TSNumber == _previousTSNumber)
		PetscFunctionReturn(0);

	// Set the previous TS number
	_previousTSNumber = TSNumber;

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

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the network
	using NetworkType = core::network::IPSIReactionNetwork;
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	// Get the number of species
	auto numSpecies = network.getSpeciesListSize();
	auto specIdI = network.getInterstitialSpeciesId();

	// Get the physical grid and step size
	auto grid = _solverHandler->getXGrid();
	double hy = _solverHandler->getStepSizeY();
	double hz = _solverHandler->getStepSizeZ();

	// Get the flux handler to know the flux amplitude.
	auto fluxHandler = _solverHandler->getFluxHandler();
	double heliumFluxAmplitude = fluxHandler->getFluxAmplitude();

	// Get the delta time from the previous timestep to this timestep
	double dt = time - _solverHandler->getPreviousTime();

	// Work of the moving surface first
	if (_solverHandler->moveSurface()) {
		// Write the initial surface positions
		if (procId == 0 && util::equal(time, 0.0)) {
			std::ofstream outputFile;
			outputFile.open("surface.txt", std::ios::app);
			outputFile << time << " ";

			// Loop on the possible yj
			for (auto yj = 0; yj < My; yj++) {
				for (auto zk = 0; zk < Mz; zk++) {
					// Get the position of the surface at yj, zk
					auto surfacePos =
						_solverHandler->getSurfacePosition(yj, zk);
					outputFile << (double)yj * hy << " " << (double)zk * hz
							   << " " << grid[surfacePos + 1] - grid[1] << " ";
				}
			}
			outputFile << std::endl;
			outputFile.close();
		}

		// Get the initial vacancy concentration
		double initialVConc = _solverHandler->getInitialVConc();

		// Loop on the possible zk and yj
		for (auto zk = 0; zk < Mz; zk++) {
			for (auto yj = 0; yj < My; yj++) {
				if (TSNumber > 0) {
					// Compute the total density of intersitials that escaped
					// from the surface since last timestep using the stored
					// flux
					_nSurf[specIdI()][yj][zk] +=
						_previousSurfFlux[specIdI()][yj][zk] * dt;

					// Remove the sputtering yield since last timestep
					_nSurf[specIdI()][yj][zk] -=
						_sputteringYield * heliumFluxAmplitude * dt;
				}

				// Get the position of the surface at yj
				auto surfacePos = _solverHandler->getSurfacePosition(yj, zk);
				auto xi = surfacePos + _solverHandler->getLeftOffset();

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
						_iClusterIds, myFlux, xi - xs);
				}

				// Gather newFlux values at this position
				double newTotalFlux = 0.0;
				MPI_Allreduce(&myFlux[specIdI()], &newTotalFlux, 1, MPI_DOUBLE,
					MPI_SUM, xolotlComm);

				// Update the previous flux
				_previousSurfFlux[specIdI()][yj][zk] = newTotalFlux;

				// Compare nInterstitials to the threshold to know if we should
				// move the surface

				// The density of tungsten is 62.8 atoms/nm3, thus the threshold
				// is
				double threshold =
					(62.8 - initialVConc) * (grid[xi] - grid[xi - 1]) * hy * hz;
				if (_nSurf[specIdI()][yj][zk] > threshold) {
					// The surface is moving
					fvalue[0] = 0.0;
				}

				// Moving the surface back
				else if (_nSurf[specIdI()][yj][zk] < -threshold / 10.0) {
					// The surface is moving
					fvalue[0] = 0.0;
				}
			}
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor3D::postEventFunction(TS ts, PetscInt nevents, PetscInt eventList[],
	PetscReal time, Vec solution, PetscBool)
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

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the network
	using NetworkType = core::network::IPSIReactionNetwork;
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	auto dof = network.getDOF();
	// Get the number of species
	auto numSpecies = network.getSpeciesListSize();
	auto specIdI = network.getInterstitialSpeciesId();

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();
	// Get the step sizes
	double hy = _solverHandler->getStepSizeY();
	double hz = _solverHandler->getStepSizeZ();

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
	double initialVConc = _solverHandler->getInitialVConc();

	// Loop on the possible zk and yj
	for (auto zk = 0; zk < Mz; zk++) {
		for (auto yj = 0; yj < My; yj++) {
			// Get the position of the surface at yj
			auto surfacePos = _solverHandler->getSurfacePosition(yj, zk);
			auto xi = surfacePos + _solverHandler->getLeftOffset();

			// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
			double threshold =
				(62.8 - initialVConc) * (grid[xi] - grid[xi - 1]) * hy * hz;

			// Move the surface up
			if (_nSurf[specIdI()][yj][zk] > threshold) {
				int nGridPoints = 0;
				// Move the surface up until it is smaller than the next
				// threshold
				while (_nSurf[specIdI()][yj][zk] > threshold &&
					surfacePos + _solverHandler->getLeftOffset() - 2 >= 0) {
					// Move the surface higher
					surfacePos--;
					xi = surfacePos + _solverHandler->getLeftOffset();
					nGridPoints++;
					// Update the number of interstitials
					_nSurf[specIdI()][yj][zk] -= threshold;
					// Update the thresold
					threshold = (62.8 - initialVConc) *
						(grid[xi] - grid[xi - 1]) * hy * hz;
				}

				// Throw an exception if the position is negative
				if (surfacePos + _solverHandler->getLeftOffset() < 2) {
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
					XOLOTL_LOG << "Adding " << nGridPoints
							   << " points to the grid on " << yj * hy
							   << " and " << zk * hz << " at time: " << time
							   << " s.";
				}

				// Set it in the solver
				_solverHandler->setSurfacePosition(surfacePos, yj, zk);

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
			else if (_nSurf[specIdI()][yj][zk] < -threshold / 10.0) {
				// Move it back as long as the number of interstitials in
				// negative
				while (_nSurf[specIdI()][yj][zk] < 0.0) {
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
					xi = surfacePos + _solverHandler->getLeftOffset();
					// Update the number of interstitials
					_nSurf[specIdI()][yj][zk] += threshold;
				}

				// Printing information about the extension of the material
				if (procId == 0) {
					XOLOTL_LOG << "Removing grid points to the grid on "
							   << yj * hy << " and " << zk * hz
							   << " at time: " << time << " s.";
				}

				// Set it in the solver
				_solverHandler->setSurfacePosition(surfacePos, yj, zk);
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
				auto surfacePos = _solverHandler->getSurfacePosition(yj, zk);
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

PetscErrorCode
PetscMonitor3D::monitorSurfaceXY(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
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

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the physical grid in the x direction
	auto grid = _solverHandler->getXGrid();

	// Setup step size variables
	double hy = _solverHandler->getStepSizeY();

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
		_surfacePlotXY->getDataProvider()->setDataPoints(myPoints);

		// Change the title of the plot and the name of the data
		std::stringstream title;
		title << "First Cluster";
		_surfacePlotXY->getDataProvider()->setDataName(title.str());
		title << " concentration";
		_surfacePlotXY->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		_surfacePlotXY->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
					  << "s";
		_surfacePlotXY->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "surfaceXY_TS" << timestep << ".png";
		_surfacePlotXY->write(fileName.str());
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor3D::monitorSurfaceXZ(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
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

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the physical grid in the x direction
	auto grid = _solverHandler->getXGrid();

	// Setup step size variables
	double hz = _solverHandler->getStepSizeZ();

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
		_surfacePlotXZ->getDataProvider()->setDataPoints(myPoints);

		// Change the title of the plot and the name of the data
		std::stringstream title;
		title << "First Cluster";
		_surfacePlotXZ->getDataProvider()->setDataName(title.str());
		title << " concentration";
		_surfacePlotXZ->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		_surfacePlotXZ->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
					  << "s";
		_surfacePlotXZ->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "surfaceXZ_TS" << timestep << ".png";
		_surfacePlotXZ->write(fileName.str());
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
} // namespace monitor
} // namespace solver
} // namespace xolotl
