#include <xolotl/core/Constants.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/io/XFile.h>
#include <xolotl/perf/ScopedTimer.h>
#include <xolotl/solver/PetscSolver.h>
#include <xolotl/solver/monitor/PetscMonitor2D.h>
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
monitorSurface(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<PetscMonitor2D*>(ictx)->monitorSurface(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

void
PetscMonitor2D::setup()
{
	PetscErrorCode ierr;

	_gbTimer = _solverHandler->getPerfHandler()->getTimer("monitor2D:GB");

	// Get the process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the xolotlViz handler registry
	auto vizHandlerRegistry = _solverHandler->getVizHandler();

	// Flags to launch the monitors or not
	PetscBool flagCheck, flagPerf, flagHeRetention, flagXeRetention, flagStatus,
		flag2DPlot, flagLargest;

	// Check the option -check_collapse
	ierr = PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCheck);
	checkPetscError(ierr,
		"setupPetsc2DMonitor: PetscOptionsHasName (-check_collapse) failed.");

	// Check the option -plot_perf
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_perf", &flagPerf);
	checkPetscError(
		ierr, "setupPetsc2DMonitor: PetscOptionsHasName (-plot_perf) failed.");

	// Check the option -plot_2d
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_2d", &flag2DPlot);
	checkPetscError(
		ierr, "setupPetsc2DMonitor: PetscOptionsHasName (-plot_2d) failed.");

	// Check the option -helium_retention
	ierr =
		PetscOptionsHasName(NULL, NULL, "-helium_retention", &flagHeRetention);
	checkPetscError(ierr,
		"setupPetsc2DMonitor: PetscOptionsHasName (-helium_retention) failed.");

	// Check the option -xenon_retention
	ierr =
		PetscOptionsHasName(NULL, NULL, "-xenon_retention", &flagXeRetention);
	checkPetscError(ierr,
		"setupPetsc2DMonitor: PetscOptionsHasName (-xenon_retention) failed.");

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, NULL, "-start_stop", &flagStatus);
	checkPetscError(
		ierr, "setupPetsc2DMonitor: PetscOptionsHasName (-start_stop) failed.");

	// Check the option -largest_conc
	ierr = PetscOptionsHasName(NULL, NULL, "-largest_conc", &flagLargest);
	checkPetscError(ierr,
		"setupPetsc2DMonitor: PetscOptionsHasName (-largest_conc) failed.");

	// Get the da from _ts
	DM da;
	ierr = TSGetDM(_ts, &da);
	checkPetscError(ierr, "setupPetsc2DMonitor: TSGetDM failed.");

	// Get the total size of the grid
	PetscInt Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
	checkPetscError(ierr, "setupPetsc2DMonitor: DMDAGetInfo failed.");

	// Get the network and its size
	auto& network = _solverHandler->getNetwork();
	// Get the number of species
	auto numSpecies = network.getSpeciesListSize();

	// Create data depending on the boundary conditions
	if (_solverHandler->getLeftOffset() == 1) {
		_nSurf = std::vector<std::vector<double>>(
			numSpecies, std::vector<double>(My, 0.0));
		_previousSurfFlux = std::vector<std::vector<double>>(
			numSpecies, std::vector<double>(My, 0.0));
	}
	if (_solverHandler->getRightOffset() == 1) {
		_nBulk = std::vector<std::vector<double>>(
			numSpecies, std::vector<double>(My, 0.0));
		_previousBulkFlux = std::vector<std::vector<double>>(
			numSpecies, std::vector<double>(My, 0.0));
	}

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
			"setupPetsc2DMonitor: PetscOptionsGetReal (-check_collapse) "
			"failed.");
		if (!flag)
			timeStepThreshold = 1.0e-16;

		// Set the post step process that tells the solver when to stop if the
		// time step collapse
		ierr = TSSetPostStep(_ts, checkTimeStep);
		checkPetscError(
			ierr, "setupPetsc2DMonitor: TSSetPostStep (checkTimeStep) failed.");
	}

	// If the user wants the surface to be able to move or bursting
	if (_solverHandler->moveSurface() || _solverHandler->burstBubbles()) {
		// Surface
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
					_nSurf[i] = lastTsGroup->readData2D(nName.str());

					// Create the previous flux attribute name
					std::ostringstream prevFluxName;
					prevFluxName << "previousFlux" << names[i] << "Surf";
					// Read the attribute
					_previousSurfFlux[i] =
						lastTsGroup->readData2D(prevFluxName.str());
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
		}

		// Bursting
		if (_solverHandler->burstBubbles()) {
			// No need to seed the random number generator here.
			// The solver handler has already done it.
		}

		// Set directions and terminate flags for the surface event
		PetscInt direction[2];
		PetscBool terminate[2];
		direction[0] = 0, direction[1] = 0;
		terminate[0] = PETSC_FALSE, terminate[1] = PETSC_FALSE;
		// Set the TSEvent
		ierr = TSSetEventHandler(_ts, 2, direction, terminate,
			monitor::eventFunction, monitor::postEventFunction, this);
		checkPetscError(ierr,
			"setupPetsc2DMonitor: TSSetEventHandler (eventFunction2D) failed.");
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
			ierr, "setupPetsc2DMonitor: TSMonitorSet (monitorPerf) failed.");
	}

	// Set the monitor to compute the helium retention
	if (flagHeRetention) {
		auto fluxHandler = _solverHandler->getFluxHandler();
		// Get the previous time if concentrations were stored and initialize
		// the fluence
		if (hasConcentrations) {
			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			_solverHandler->setPreviousTime(previousTime);
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
					_nSurf[i] = lastTsGroup->readData2D(nName.str());

					// Create the previous flux attribute name
					std::ostringstream prevFluxName;
					prevFluxName << "previousFlux" << names[i] << "Surf";
					// Read the attribute
					_previousSurfFlux[i] =
						lastTsGroup->readData2D(prevFluxName.str());
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
					_nBulk[i] = lastTsGroup->readData2D(nName.str());

					// Create the previous flux attribute name
					std::ostringstream prevFluxName;
					prevFluxName << "previousFlux" << names[i] << "Bulk";
					// Read the attribute
					_previousBulkFlux[i] =
						lastTsGroup->readData2D(prevFluxName.str());
				}
			}

			// Bursting
			if (_solverHandler->burstBubbles()) {
				// Read about the impurity fluxes in from bursting
				_nHeliumBurst = lastTsGroup->readData1D("nHeliumBurst");
				_nDeuteriumBurst = lastTsGroup->readData1D("nDeuteriumBurst");
				_nTritiumBurst = lastTsGroup->readData1D("nTritiumBurst");
			}
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeFluence, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc2DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeHeliumRetention2D will be called at each timestep
		ierr =
			TSMonitorSet(_ts, monitor::computeHeliumRetention, this, nullptr);
		checkPetscError(ierr,
			"setupPetsc2DMonitor: TSMonitorSet (computeHeliumRetention) "
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
			outputFile << "Helium_burst Deuterium_burst Tritium_burst"
					   << std::endl;
			outputFile.close();

			if (_solverHandler->temporalFlux()) {
				// Open an additional file that will keep the flux evolution
				outputFile.open("instantFlux.txt");
				outputFile << "#time ";

				// Get the generated clusters
				auto indices = fluxHandler->getFluxIndices();

				// Get the bounds
				auto bounds = network.getAllClusterBounds();
				// Loop on them
				for (auto i : indices) {
					for (auto id = core::network::SpeciesId(numSpecies); id;
						 ++id) {
						auto speciesName = network.getSpeciesName(id);
						if (bounds[i][2 * id()] > 0)
							outputFile << speciesName << "_"
									   << bounds[i][2 * id()];
					}
					outputFile << " ";
				}
				outputFile << std::endl;
				outputFile.close();
			}
		}
	}

	// Set the monitor to compute the xenon retention
	if (flagXeRetention) {
		// Get the da from _ts
		DM da;
		ierr = TSGetDM(_ts, &da);
		checkPetscError(ierr, "setupPetsc2DMonitor: TSGetDM failed.");
		// Get the local boundaries
		PetscInt xm, ym;
		ierr = DMDAGetCorners(da, NULL, NULL, NULL, &xm, &ym, NULL);
		checkPetscError(ierr, "setupPetsc2DMonitor: DMDAGetCorners failed.");
		// Create the local vectors on each process
		_solverHandler->createLocalNE(xm, ym);

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
			ierr, "setupPetsc2DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeXenonRetention2D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeXenonRetention, this, nullptr);
		checkPetscError(ierr,
			"setupPetsc2DMonitor: TSMonitorSet (computeXenonRetention) "
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
	if (flag2DPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a SurfacePlot
			_surfacePlot = vizHandlerRegistry->getPlot(viz::PlotType::SURFACE);

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>();
			labelProvider->axis1Label = "Depth (nm)";
			labelProvider->axis2Label = "Y (nm)";
			labelProvider->axis3Label = "Concentration";

			// Give it to the plot
			_surfacePlot->setLabelProvider(labelProvider);

			// Create the data provider
			auto dataProvider =
				std::make_shared<viz::dataprovider::CvsXYDataProvider>();

			// Give it to the plot
			_surfacePlot->setDataProvider(dataProvider);
		}

		// monitorSurface2D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorSurface, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc2DMonitor: TSMonitorSet (monitorSurface) failed.");
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
			"setupPetsc2DMonitor: PetscOptionsGetReal (-largest_conc) failed.");

		// monitorLargest2D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorLargest, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc2DMonitor: TSMonitorSet (monitorLargest) failed.");
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr =
			PetscOptionsGetReal(NULL, NULL, "-start_stop", &_hdf5Stride, &flag);
		checkPetscError(ierr,
			"setupPetsc2DMonitor: PetscOptionsGetReal (-start_stop) failed.");
		if (!flag)
			_hdf5Stride = 1.0;

		if (hasConcentrations) {
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

			// Create and initialize a checkpoint file.
			// We do this in its own scope so that the file
			// is closed when the file object goes out of scope.
			// We want it to close before we (potentially) copy
			// the network from another file using a single-process
			// MPI communicator.
			{
				io::XFile checkpointFile(
					_hdf5OutputName, grid, xolotlComm, My, hy);
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

		// startStop2D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::startStop, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc2DMonitor: TSMonitorSet (startStop) failed.");
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(_ts, monitor::monitorTime, this, nullptr);
	checkPetscError(
		ierr, "setupPetsc2DMonitor: TSMonitorSet (monitorTime) failed.");
}

PetscErrorCode
PetscMonitor2D::monitorLargest(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declaration
	PetscErrorCode ierr;
	double ***solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Loop on the local grid
	for (auto j = ys; j < ys + ym; j++)
		for (auto i = xs; i < xs + xm; i++) {
			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[j][i];
			// Check the concentration
			if (gridPointSolution[_largestClusterId] > _largestThreshold) {
				ierr = TSSetConvergedReason(ts, TS_CONVERGED_USER);
				CHKERRQ(ierr);
				// Send an error
				throw std::runtime_error(
					"\nxolotlSolver::Monitor2D: The largest "
					"cluster concentration is too high!!");
			}
		}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor2D::startStop(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declaration
	PetscErrorCode ierr;
	const double ***solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Get local coordinates
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
	std::vector<int> surfaceIndices;
	for (auto i = 0; i < My; i++) {
		surfaceIndices.push_back(_solverHandler->getSurfacePosition(i));
	}

	// Open the existing checkpoint file.
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
		tsGroup->writeSurface2D(
			surfaceIndices, _nSurf, _previousSurfFlux, names);
	}

	// Write the bottom impurity information if the bottom is a free surface
	if (_solverHandler->getRightOffset() == 1)
		tsGroup->writeBottom2D(_nBulk, _previousBulkFlux, names);

	// Write the bursting information if the bubble bursting is used
	if (_solverHandler->burstBubbles())
		tsGroup->writeBursting(_nHeliumBurst, _nDeuteriumBurst, _nTritiumBurst);

	// Loop on the full grid
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
			if (i >= xs && i < xs + xm && j >= ys && j < ys + ym) {
				write = true;
				// Get the pointer to the beginning of the solution data for
				// this grid point
				gridPointSolution = solutionArray[j][i];

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
			MPI_Allreduce(&concId, &concProc, 1, MPI_INT, MPI_SUM, xolotlComm);

			// Broadcast the size
			MPI_Bcast(&concSize, 1, MPI_INT, concProc, xolotlComm);

			// Skip the grid point if the size is 0
			if (concSize == 0)
				continue;

			// All processes create the dataset and fill it
			tsGroup->writeConcentrationDataset(
				concSize, concArray, write, i, j);
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor2D::computeHeliumRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the flux handler.
	auto fluxHandler = _solverHandler->getFluxHandler();
	// Get the diffusion handler
	auto diffusionHandler = _solverHandler->getDiffusionHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the physical grid in the x direction
	auto grid = _solverHandler->getXGrid();

	// Setup step size variables
	double hy = _solverHandler->getStepSizeY();

	// Get the network
	using NetworkType = core::network::IPSIReactionNetwork;
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	const auto dof = network.getDOF();

	// Get the array of concentration
	double ***solutionArray, *gridPointSolution;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration over the grid
	auto numSpecies = network.getSpeciesListSize();
	auto specIdI = network.getInterstitialSpeciesId();
	auto myConcData = std::vector<double>(numSpecies, 0.0);

	// Loop on the grid
	for (auto yj = ys; yj < ys + ym; yj++) {
		// Get the surface position
		auto surfacePos = _solverHandler->getSurfacePosition(yj);

		for (auto xi = xs; xi < xs + xm; xi++) {
			// Boundary conditions
			if (xi < surfacePos + _solverHandler->getLeftOffset() ||
				xi >= Mx - _solverHandler->getRightOffset())
				continue;

			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[yj][xi];

			double hx = grid[xi + 1] - grid[xi];

			using HostUnmanaged = Kokkos::View<double*, Kokkos::HostSpace,
				Kokkos::MemoryUnmanaged>;
			auto hConcs = HostUnmanaged(gridPointSolution, dof);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof);
			deep_copy(dConcs, hConcs);

			// Get the total concentrations at this grid point
			for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
				myConcData[id()] +=
					network.getTotalAtomConcentration(dConcs, id, 1) * hx * hy;
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

		// Loop on every Y position
		for (auto j = 0; j < My; j++) {
			// Set the surface position
			auto xi = _solverHandler->getSurfacePosition(j) + 1;
			// Value to know on which processor is the bottom
			int surfProc = 0;

			// Check we are on the right proc
			if (xi >= xs && xi < xs + xm && j >= ys && j < ys + ym) {
				// Compute the total number of impurities that left at the
				// surface
				if (timestep > 0) {
					for (auto i = 0; i < numSpecies; ++i) {
						if (i == specIdI() && _solverHandler->moveSurface())
							continue;
						_nSurf[i][j] += _previousSurfFlux[i][j] * dt;
					}
				}
				auto myFluxData = std::vector<double>(numSpecies, 0.0);

				// Get the pointer to the beginning of the solution data for
				// this grid point
				gridPointSolution = solutionArray[j][xi];

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
				double factor = 2.0 * hy / (hxLeft + hxRight);

				network.updateOutgoingDiffFluxes(gridPointSolution, factor,
					diffusingIds, myFluxData, xi - xs);

				// Take into account the surface advection
				// Get the surface advection handler
				auto advecHandler = _solverHandler->getAdvectionHandler();
				// Get the sink strengths and advecting clusters
				auto sinkStrengths = advecHandler->getSinkStrengths();
				auto advecClusters = advecHandler->getAdvectingClusters();
				// Set the distance from the surface
				double distance = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1] -
					advecHandler->getLocation();

				network.updateOutgoingAdvecFluxes(gridPointSolution,
					3.0 * hy /
						(core::kBoltzmann * distance * distance * distance *
							distance),
					advecClusters, sinkStrengths, myFluxData, xi - xs);

				for (auto i = 0; i < numSpecies; ++i) {
					if (i == specIdI() && _solverHandler->moveSurface())
						continue;
					_previousSurfFlux[i][j] = myFluxData[i];
				}

				// Set the bottom processor
				surfProc = procId;
			}

			// Get which processor will send the information
			int surfId = 0;
			MPI_Allreduce(&surfProc, &surfId, 1, MPI_INT, MPI_SUM, xolotlComm);

			// Send the information about impurities
			// to the other processes
			std::vector<double> countFluxData;
			for (auto i = 0; i < numSpecies; ++i) {
				countFluxData.push_back(_nSurf[i][j]);
				countFluxData.push_back(_previousSurfFlux[i][j]);
			}
			MPI_Bcast(countFluxData.data(), countFluxData.size(), MPI_DOUBLE,
				surfId, xolotlComm);

			// Extract impurity data from broadcast buffer.
			for (auto i = 0; i < numSpecies; ++i) {
				_nSurf[i][j] = countFluxData[2 * i];
				_previousSurfFlux[i][j] = countFluxData[(2 * i) + 1];
			}
		}
	}

	// Look at the fluxes going in the bulk if the bottom is a free surface
	if (_solverHandler->getRightOffset() == 1) {
		// Set the bottom surface position
		auto xi = Mx - 2;

		// Get the vector of diffusing clusters
		auto diffusingIds = diffusionHandler->getDiffusingIds();

		// Loop on every Y position
		for (auto j = 0; j < My; j++) {
			// Value to know on which processor is the bottom
			int bottomProc = 0;

			// Check we are on the right proc
			if (xi >= xs && xi < xs + xm && j >= ys && j < ys + ym) {
				// Compute the total number of impurities that left at the
				// surface
				if (timestep > 0) {
					for (auto i = 0; i < numSpecies; ++i) {
						_nBulk[i][j] += _previousBulkFlux[i][j] * dt;
					}
				}
				auto myFluxData = std::vector<double>(numSpecies, 0.0);

				// Get the pointer to the beginning of the solution data for
				// this grid point
				gridPointSolution = solutionArray[j][xi];

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
				double factor = 2.0 * hy / (hxLeft + hxRight);

				network.updateOutgoingDiffFluxes(gridPointSolution, factor,
					diffusingIds, myFluxData, xi - xs);

				for (auto i = 0; i < numSpecies; ++i) {
					_previousBulkFlux[i][j] = myFluxData[i];
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
				countFluxData.push_back(_nBulk[i][j]);
				countFluxData.push_back(_previousBulkFlux[i][j]);
			}
			MPI_Bcast(countFluxData.data(), countFluxData.size(), MPI_DOUBLE,
				bottomId, xolotlComm);

			// Extract impurity data from broadcast buffer.
			for (auto i = 0; i < numSpecies; ++i) {
				_nBulk[i][j] = countFluxData[2 * i];
				_previousBulkFlux[i][j] = countFluxData[(2 * i) + 1];
			}
		}
	}

	// Master process
	if (procId == 0) {
		// Compute the total surface irradiated by the helium flux
		double surface = (double)My * hy;

		// Rescale the concentration
		for (auto i = 0; i < numSpecies; ++i) {
			totalConcData[i] /= surface;
		}
		auto totalBulk = std::vector<double>(numSpecies, 0.0);
		// Look if the bottom is a free surface
		if (_solverHandler->getRightOffset() == 1) {
			for (auto i = 0; i < numSpecies; ++i) {
				for (auto j = 0; j < My; j++) {
					totalBulk[i] += _nBulk[i][j];
				}
				totalBulk[i] = totalBulk[i] / surface;
			}
		}
		auto totalSurf = std::vector<double>(numSpecies, 0.0);
		// Look if the bottom is a free surface
		if (_solverHandler->getLeftOffset() == 1) {
			for (auto i = 0; i < numSpecies; ++i) {
				for (auto j = 0; j < My; j++) {
					totalSurf[i] += _nSurf[i][j];
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
		outputFile << _nHeliumBurst << " " << _nDeuteriumBurst << " "
				   << _nTritiumBurst << std::endl;
		outputFile.close();

		if (_solverHandler->temporalFlux()) {
			// Open an additional file that will keep the flux evolution
			outputFile.open("instantFlux.txt", std::ios::app);
			outputFile << time << " ";
			// Get the flux information
			auto instantFlux = fluxHandler->getInstantFlux(time);
			// Loop on it
			for (auto flux : instantFlux) {
				outputFile << flux << " ";
			}

			outputFile << std::endl;
			outputFile.close();
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor2D::computeXenonRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	perf::ScopedTimer myTimer(_gbTimer);

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
	PetscReal ***solutionArray, *gridPointSolution;
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
	auto xeCluster = network.findCluster(xeComp, plsm::HostMemSpace{});
	auto xeId = xeCluster.getId();

	// Loop on the grid
	for (auto yj = ys; yj < ys + ym; yj++)
		for (auto xi = xs; xi < xs + xm; xi++) {
			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[yj][xi];

			using HostUnmanaged = Kokkos::View<double*, Kokkos::HostSpace,
				Kokkos::MemoryUnmanaged>;
			auto hConcs = HostUnmanaged(gridPointSolution, dof);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof);
			deep_copy(dConcs, hConcs);

			double hx = grid[xi + 1] - grid[xi];

			// Get the concentrations
			xeConcentration +=
				network.getTotalAtomConcentration(dConcs, Spec::Xe, 1) * hx *
				hy;
			bubbleConcentration +=
				network.getTotalConcentration(dConcs, Spec::Xe, 1) * hx * hy;
			radii += network.getTotalRadiusConcentration(dConcs, Spec::Xe, 1) *
				hx * hy;
			partialBubbleConcentration +=
				network.getTotalConcentration(dConcs, Spec::Xe, minSizes[0]) *
				hx * hy;
			partialRadii += network.getTotalRadiusConcentration(
								dConcs, Spec::Xe, minSizes[0]) *
				hx * hy;

			// Set the volume fraction
			double volumeFrac =
				network.getTotalVolumeFraction(dConcs, Spec::Xe, minSizes[0]);
			_solverHandler->setVolumeFraction(volumeFrac, xi - xs, yj - ys);
			// Set the monomer concentration
			_solverHandler->setMonomerConc(
				gridPointSolution[xeCluster.getId()], xi - xs, yj - ys);
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
	// Get the previous Xe flux vector
	auto& localNE = _solverHandler->getLocalNE();
	// Loop on the GB
	for (auto const& pair : gbVector) {
		// Middle
		auto xi = std::get<0>(pair);
		auto yj = std::get<1>(pair);
		// Check we are on the right proc
		if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
			double previousXeFlux = std::get<1>(localNE[xi - xs][yj - ys][0]);
			globalXeFlux += previousXeFlux * (grid[xi + 1] - grid[xi]) * hy;
			// Set the amount in the vector we keep
			_solverHandler->setLocalXeRate(
				previousXeFlux * dt, xi - xs, yj - ys);
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
		if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
			// X segment
			// Left
			xi = std::get<0>(pair) - 1;
			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[yj][xi];
			// Compute the flux coming from the left
			network.updateOutgoingDiffFluxes(gridPointSolution, factor / hxLeft,
				diffusingIds, myRate, xi + 1 - xs);

			// Right
			xi = std::get<0>(pair) + 1;
			gridPointSolution = solutionArray[yj][xi];
			// Compute the flux coming from the right
			network.updateOutgoingDiffFluxes(gridPointSolution,
				factor / hxRight, diffusingIds, myRate, xi + 1 - xs);

			// Y segment
			// Bottom
			xi = std::get<0>(pair);
			yj = std::get<1>(pair) - 1;
			gridPointSolution = solutionArray[yj][xi];
			// Compute the flux coming from the bottom
			network.updateOutgoingDiffFluxes(gridPointSolution, 1.0 / (hy * hy),
				diffusingIds, myRate, xi + 1 - xs);

			// Top
			yj = std::get<1>(pair) + 1;
			gridPointSolution = solutionArray[yj][xi];
			// Compute the flux coming from the top
			network.updateOutgoingDiffFluxes(gridPointSolution, 1.0 / (hy * hy),
				diffusingIds, myRate, xi + 1 - xs);

			// Middle
			xi = std::get<0>(pair);
			yj = std::get<1>(pair);
			_solverHandler->setPreviousXeFlux(myRate[0], xi - xs, yj - ys);
		}
	}

	// Master process
	if (procId == 0) {
		// Compute the total surface irradiated
		double surface = (double)My * hy;
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
PetscMonitor2D::eventFunction(
	TS ts, PetscReal time, Vec solution, PetscScalar* fvalue)
{
	// Initial declaration
	PetscErrorCode ierr;
	double ***solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	fvalue[0] = 1.0, fvalue[1] = 1.0;
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

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();
	// Get the step size in Y
	double hy = _solverHandler->getStepSizeY();

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
				// Get the position of the surface at yj
				auto surfacePos = _solverHandler->getSurfacePosition(yj);
				outputFile << grid[surfacePos + 1] - grid[1] << " ";
			}
			outputFile << std::endl;
			outputFile.close();
		}

		// Loop on the possible yj
		for (auto yj = 0; yj < My; yj++) {
			if (TSNumber > 0) {
				// Compute the total density of intersitials that escaped from
				// the surface since last timestep using the stored flux
				_nSurf[specIdI()][yj] += _previousSurfFlux[specIdI()][yj] * dt;

				// Remove the sputtering yield since last timestep
				_nSurf[specIdI()][yj] -=
					_sputteringYield * heliumFluxAmplitude * dt * hy;
			}

			// Get the position of the surface at yj
			const auto surfacePos = _solverHandler->getSurfacePosition(yj);
			auto xi = surfacePos + _solverHandler->getLeftOffset();
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

			// Initialize the value for the flux
			auto myFlux = std::vector<double>(numSpecies, 0.0);

			// if xi is on this process
			if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
				// Get the concentrations at xi = surfacePos + 1
				gridPointSolution = solutionArray[yj][xi];

				// Factor for finite difference
				double factor = hy * 2.0 / (hxLeft + hxRight);

				network.updateOutgoingDiffFluxes(
					gridPointSolution, factor, _iClusterIds, myFlux, xi - xs);
			}

			// Check if the surface on the left and/or right sides are higher
			// than at this grid point
			PetscInt yLeft = yj - 1, yRight = yj + 1;
			if (yLeft < 0)
				yLeft = My - 1; // Periodicity
			if (yRight == My)
				yRight = 0; // Periodicity
			// We want the position just at the surface now
			xi = surfacePos;
			// Do the left side first
			auto leftSurf = _solverHandler->getSurfacePosition(yLeft);
			if (leftSurf < surfacePos) {
				// Loop on every grid point above
				for (auto currentX = xi; currentX >= xi; --currentX) {
					// If this is the locally owned part of the grid
					if (currentX >= xs && currentX < xs + xm && yLeft >= ys &&
						yLeft < ys + ym) {
						// Get the concentrations at xi = surfacePos
						gridPointSolution = solutionArray[yLeft][currentX];
						double hX = 0.0;
						if (currentX - 1 < 0) {
							hX = grid[currentX + 1] - grid[currentX];
						}
						else {
							hX =
								(grid[currentX + 1] - grid[currentX - 1]) / 2.0;
						}

						network.updateOutgoingDiffFluxes(gridPointSolution,
							hX / hy, _iClusterIds, myFlux, currentX - xs);
					}
				}
			}
			// Now do the right side
			auto rightSurf = _solverHandler->getSurfacePosition(yRight);
			if (rightSurf < surfacePos) {
				// Loop on every grid point above
				for (auto currentX = xi; currentX >= xi; --currentX) {
					// If this is the locally owned part of the grid
					if (currentX >= xs && currentX < xs + xm && yRight >= ys &&
						yRight < ys + ym) {
						// Get the concentrations at xi = surfacePos
						gridPointSolution = solutionArray[yRight][currentX];
						double hX = 0.0;
						if (currentX - 1 < 0) {
							hX = grid[currentX + 1] - grid[currentX];
						}
						else {
							hX =
								(grid[currentX + 1] - grid[currentX - 1]) / 2.0;
						}

						network.updateOutgoingDiffFluxes(gridPointSolution,
							hX / hy, _iClusterIds, myFlux, currentX - xs);
					}
				}
			}

			// Gather newFlux values at this position
			double newTotalFlux = 0.0;
			MPI_Allreduce(&myFlux[specIdI()], &newTotalFlux, 1, MPI_DOUBLE,
				MPI_SUM, xolotlComm);

			// Update the previous flux
			_previousSurfFlux[specIdI()][yj] = newTotalFlux;

			// Compare nInterstitials to the threshold to know if we should move
			// the surface

			// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
			double threshold =
				core::tungstenDensity * (grid[xi] - grid[xi - 1]) * hy;
			if (_nSurf[specIdI()][yj] > threshold) {
				// The surface is moving
				fvalue[0] = 0.0;
			}

			// Update the threshold for erosion (the cell size is not the same)
			threshold = core::tungstenDensity * (grid[xi + 1] - grid[xi]) * hy;
			// Moving the surface back
			if (_nSurf[specIdI()][yj] < -threshold * 0.9) {
				// The surface is moving
				fvalue[0] = 0.0;
			}
		}
	}

	// Now work on the bubble bursting
	if (_solverHandler->burstBubbles()) {
		auto dof = network.getDOF();
		auto specIdHe = network.getHeliumSpeciesId();

		// Compute the prefactor for the probability (arbitrary)
		double prefactor =
			heliumFluxAmplitude * dt * _solverHandler->getBurstingFactor();

		// The depth parameter to know where the bursting should happen
		double depthParam = _solverHandler->getTauBursting(); // nm
		// The number of He per V in a bubble
		double heVRatio = _solverHandler->getHeVRatio();

		// For now we are not bursting
		bool burst = false;

		// Loop on the full grid
		for (auto yj = 0; yj < My; yj++) {
			// Get the surface position
			auto surfacePos = _solverHandler->getSurfacePosition(yj);
			for (auto xi = surfacePos + _solverHandler->getLeftOffset();
				 xi < Mx - _solverHandler->getRightOffset(); xi++) {
				bool localBurst = false;
				// If this is the locally owned part of the grid
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
					// Get the pointer to the beginning of the solution data for
					// this grid point
					gridPointSolution = solutionArray[yj][xi];

					using HostUnmanaged = Kokkos::View<double*,
						Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
					auto hConcs = HostUnmanaged(gridPointSolution, dof);
					auto dConcs = Kokkos::View<double*>("Concentrations", dof);
					deep_copy(dConcs, hConcs);

					// Get the distance from the surface
					double distance =
						(grid[xi] + grid[xi + 1]) / 2.0 - grid[surfacePos + 1];

					// Compute the helium density at this grid point
					double heDensity =
						network.getTotalAtomConcentration(dConcs, specIdHe, 1);

					// Compute the radius of the bubble from the number of
					// helium
					double nV =
						heDensity * (grid[xi + 1] - grid[xi]) * hy / heVRatio;
					double latticeParam = network.getLatticeParameter();
					double tlcCubed =
						latticeParam * latticeParam * latticeParam;
					double radius = (sqrt(3.0) / 4) * latticeParam +
						cbrt((3.0 * tlcCubed * nV) / (8.0 * core::pi)) -
						cbrt((3.0 * tlcCubed) / (8.0 * core::pi));

					// Add randomness
					double prob = prefactor *
						(1.0 - (distance - radius) / distance) *
						std::min(1.0,
							exp(-(distance - depthParam) / (depthParam * 2.0)));
					double test = _solverHandler->getRNG().GetRandomDouble();

					if (radius > distance || prob > test) {
						localBurst = true;
					}
				}
				// Check if this location burst
				bool locationBurst = false;
				MPI_Allreduce(&localBurst, &locationBurst, 1, MPI_C_BOOL,
					MPI_LOR, xolotlComm);
				if (locationBurst) {
					_depthPositions.push_back(std::make_pair(yj, xi));
					burst = true;
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

PetscErrorCode
PetscMonitor2D::postEventFunction(TS ts, PetscInt nevents, PetscInt eventList[],
	PetscReal time, Vec solution, PetscBool)
{
	// Initial declaration
	PetscErrorCode ierr;
	double ***solutionArray, *gridPointSolution;
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
	auto& network = _solverHandler->getNetwork();
	auto dof = network.getDOF();

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();
	// Get the step size in Y
	double hy = _solverHandler->getStepSizeY();

	// Take care of bursting
	using NetworkType = core::network::IPSIReactionNetwork;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);
	auto specIdV = psiNetwork->getVacancySpeciesId();
	auto specIdI = psiNetwork->getInterstitialSpeciesId();
	auto nBurst = std::vector<double>(3, 0.0);

	bool surfaceMoved = false;

	// Loop on each bursting depth
	for (auto i = 0; i < _depthPositions.size(); i++) {
		// Get the coordinates of the point
		auto xi = _depthPositions[i].second, yj = _depthPositions[i].first;
		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[yj][xi];

		// Get the surface position
		auto surfacePos = _solverHandler->getSurfacePosition(yj);
		// Get the distance from the surface
		double distance =
			(grid[xi] + grid[xi + 1]) / 2.0 - grid[surfacePos + 1];

		std::cout << "bursting at: " << yj * hy << " " << distance << std::endl;

		// Pinhole case
		if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
			psiNetwork->updateBurstingConcs(gridPointSolution, 0.0, nBurst);
		}

		//		// Crater case
		//		PetscInt yLeft = yj - 1, yRight = yj + 1;
		//		if (yLeft < 0)
		//			yLeft = My - 1; // Periodicity
		//		if (yRight == My)
		//			yRight = 0; // Periodicity
		//
		//		// Loop on every grid point above
		//		for (auto currentX = xi; currentX > surfacePos; --currentX) {
		//			// If this is the locally owned part of the grid
		//			if (currentX >= xs && currentX < xs + xm && yj >= ys &&
		//				yj < ys + ym) {
		//				gridPointSolution = solutionArray[yj][currentX];
		//				// Get the total I and V concentrations
		//				using HostUnmanaged = Kokkos::View<double*,
		// Kokkos::HostSpace, 					Kokkos::MemoryUnmanaged>;
		// auto hConcs
		// = HostUnmanaged(gridPointSolution, dof); 				auto dConcs
		// = Kokkos::View<double*>("Concentrations", dof); deep_copy(dConcs,
		// hConcs); 				double iConc =
		// psiNetwork->getTotalAtomConcentration(dConcs, specIdI, 1);
		// double vConc =
		// psiNetwork->getTotalAtomConcentration(dConcs, specIdV, 1);
		//				// The density of tungsten is 62.8 atoms/nm3
		//				double wConc = (core::tungstenDensity - vConc + iConc)
		/// 2.0;
		//
		//				// Reset the concentrations
		//				for (auto l = 0; l < dof; ++l) {
		//					gridPointSolution[l] = 0.0;
		//				}
		//
		//				// Pass the tungsten concentration to the sides
		//				// Do the left side first
		//				if (_solverHandler->getSurfacePosition(yLeft) <
		// currentX)
		//{
		//					// if we are on the right process
		//					if (currentX >= xs && currentX < xs + xm && yLeft >=
		// ys
		//&& 						yLeft < ys + ym) {
		//						// Get the concentrations at currentX
		//						gridPointSolution =
		// solutionArray[yLeft][currentX];
		// gridPointSolution[_iClusterIds[0]] += wConc;
		//					}
		//				}
		//				// Now do the right side
		//				if (_solverHandler->getSurfacePosition(yRight) <
		// currentX)
		//{
		//					// if we are on the right process
		//					if (currentX >= xs && currentX < xs + xm && yRight
		//>= ys
		//&& 						yRight < ys + ym) {
		//						// Get the concentrations at currentX
		//						gridPointSolution =
		// solutionArray[yRight][currentX];
		// gridPointSolution[_iClusterIds[0]]
		//+= wConc;
		//					}
		//				}
		//			}
		//		}
		//
		//		// Update the surface position
		//		_solverHandler->setSurfacePosition(xi, yj);
		//
		//		// Reset the I flux and count
		//		nInterstitial2D[yLeft] += nInterstitial2D[yj] / 2.0;
		//		nInterstitial2D[yRight] += nInterstitial2D[yj] / 2.0;
		//		nInterstitial2D[yj] = 0.0;
		//		previousIFlux2D[yj] = 0.0;
		//
		//		surfaceMoved = true;
	}

	// Add up the local quantities
	auto globalBurst = std::vector<double>(3, 0.0);
	MPI_Allreduce(
		nBurst.data(), globalBurst.data(), 3, MPI_DOUBLE, MPI_SUM, xolotlComm);
	_nHeliumBurst += globalBurst[0];
	_nDeuteriumBurst += globalBurst[1];
	_nTritiumBurst += globalBurst[2];

	// Now takes care of moving surface
	bool moving = false;
	for (auto i = 0; i < nevents; i++) {
		if (eventList[i] == 0)
			moving = true;
	}

	// Skip if nothing is moving
	if (!moving && !surfaceMoved) {
		// Restore the solutionArray
		ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
		CHKERRQ(ierr);

		PetscFunctionReturn(0);
	}

	// Loop on the possible yj
	for (auto yj = 0; yj < My; yj++) {
		// Get the position of the surface at yj
		auto surfacePos = _solverHandler->getSurfacePosition(yj);
		auto xi = surfacePos + _solverHandler->getLeftOffset();

		// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
		double threshold =
			core::tungstenDensity * (grid[xi] - grid[xi - 1]) * hy;

		// Move the surface up
		if (_nSurf[specIdI()][yj] > threshold) {
			int nGridPoints = 0;
			// Move the surface up until it is smaller than the next threshold
			while (_nSurf[specIdI()][yj] > threshold &&
				surfacePos + _solverHandler->getLeftOffset() - 2 >= 0) {
				// Move the surface higher
				surfacePos--;
				xi = surfacePos + _solverHandler->getLeftOffset();
				nGridPoints++;
				// Update the number of interstitials
				_nSurf[specIdI()][yj] -= threshold;
				// Update the thresold
				threshold =
					core::tungstenDensity * (grid[xi] - grid[xi - 1]) * hy;
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
					"\nxolotlSolver::Monitor2D: The surface is "
					"trying to go outside of the grid!!");
			}

			// Printing information about the extension of the material
			if (procId == 0) {
				XOLOTL_LOG << "Adding " << nGridPoints
						   << " points to the grid on " << yj * hy
						   << " at time: " << time << " s.";
			}

			// Set it in the solver
			_solverHandler->setSurfacePosition(surfacePos, yj);

			// Initialize the vacancy concentration and the temperature on the
			// new grid points Get the single vacancy ID
			auto singleVacancyCluster = network.getSingleVacancy();
			auto vacancyIndex = core::network::IReactionNetwork::invalidIndex();
			if (singleVacancyCluster.getId() !=
				core::network::IReactionNetwork::invalidIndex())
				vacancyIndex = singleVacancyCluster.getId();
			// Get the surface temperature
			double temp = 0.0;
			if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
				temp = solutionArray[yj][xi][dof];
			}
			double surfTemp = 0.0;
			MPI_Allreduce(&temp, &surfTemp, 1, MPI_DOUBLE, MPI_SUM, xolotlComm);
			// Loop on the new grid points
			while (nGridPoints >= 0) {
				// Position of the newly created grid point
				xi = surfacePos + nGridPoints;

				// If xi is on this process
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
					// Get the concentrations
					gridPointSolution = solutionArray[yj][xi];

					// Set the new surface temperature
					gridPointSolution[dof] = surfTemp;

					// Reset the concentrations
					for (auto l = 0; l < dof; ++l) {
						gridPointSolution[l] = 0.0;
					}

					auto initialConc = _solverHandler->getInitialConc();

					if (vacancyIndex > 0 && nGridPoints > 0) {
						// Initialize the concentration
						for (auto pair : initialConc) {
							gridPointSolution[pair.first] = pair.second;
						}
					}
				}

				// Decrease the number of grid points
				--nGridPoints;
			}
		}

		// Moving the surface back
		else if (_nSurf[specIdI()][yj] < -threshold / 10.0) {
			// Move it back as long as the number of interstitials in negative
			while (_nSurf[specIdI()][yj] < 0.0) {
				// Compute the threshold to a deeper grid point
				threshold =
					core::tungstenDensity * (grid[xi + 1] - grid[xi]) * hy;
				// Set all the concentrations to 0.0 at xi = surfacePos + 1
				// if xi is on this process
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
					// Get the concentrations at xi = surfacePos + 1
					gridPointSolution = solutionArray[yj][xi];
					// Loop on DOF
					for (auto i = 0; i < dof; i++) {
						gridPointSolution[i] = 0.0;
					}
				}

				// Move the surface deeper
				surfacePos++;
				xi = surfacePos + _solverHandler->getLeftOffset();
				// Update the number of interstitials
				_nSurf[specIdI()][yj] += threshold;
			}

			// Printing information about the extension of the material
			if (procId == 0) {
				XOLOTL_LOG << "Removing grid points to the grid on " << yj * hy
						   << " at time: " << time << " s.";
			}

			// Set it in the solver
			_solverHandler->setSurfacePosition(surfacePos, yj);
		}
	}

	// Write the surface positions
	if (procId == 0) {
		std::ofstream outputFile;
		outputFile.open("surface.txt", std::ios::app);
		outputFile << time << " ";

		// Loop on the possible yj
		for (auto yj = 0; yj < My; yj++) {
			// Get the position of the surface at yj
			auto surfacePos = _solverHandler->getSurfacePosition(yj);
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

PetscErrorCode
PetscMonitor2D::monitorSurface(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	const double ***solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	double x = 0.0, y = 0.0;

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

	// Loop on the full grid
	for (auto j = 0; j < My; j++) {
		for (auto i = 0; i < Mx; i++) {
			// If it is the locally owned part of the grid
			if (i >= xs && i < xs + xm && j >= ys && j < ys + ym) {
				// Get the pointer to the beginning of the solution data for
				// this grid point
				gridPointSolution = solutionArray[j][i];
				// Compute x and y
				x = (grid[i] + grid[i + 1]) / 2.0 - grid[1];
				y = (double)j * hy;

				// If it is procId 0 just store the value in the myPoints vector
				if (procId == 0) {
					// Modify the Point with the gridPointSolution[iCluster] as
					// the value and add it to myPoints
					thePoint.value = gridPointSolution[iCluster];
					thePoint.t = time;
					thePoint.x = x;
					thePoint.y = y;
					myPoints->push_back(thePoint);
				}
				// Else, the values must be sent to procId 0
				else {
					// Send the value of the local position to the master
					// process
					MPI_Send(&x, 1, MPI_DOUBLE, 0, 10, xolotlComm);
					// Send the value of the local position to the master
					// process
					MPI_Send(&y, 1, MPI_DOUBLE, 0, 11, xolotlComm);

					// Send the value of the concentration to the master process
					MPI_Send(&gridPointSolution[iCluster], 1, MPI_DOUBLE, 0, 12,
						xolotlComm);
				}
			}
			// Else if it is NOT the locally owned part of the grid but still
			// procId == 0, it should receive the values for the point and add
			// them to myPoint
			else if (procId == 0) {
				// Get the position
				MPI_Recv(&x, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 10, xolotlComm,
					MPI_STATUS_IGNORE);
				MPI_Recv(&y, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 11, xolotlComm,
					MPI_STATUS_IGNORE);

				// and the concentration
				double conc = 0.0;
				MPI_Recv(&conc, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 12, xolotlComm,
					MPI_STATUS_IGNORE);

				// Modify the Point with the received values and add it to
				// myPoints
				thePoint.value = conc;
				thePoint.t = time;
				thePoint.x = x;
				thePoint.y = y;
				myPoints->push_back(thePoint);
			}

			// Wait for everybody at each grid point
			MPI_Barrier(xolotlComm);
		}
	}

	// Plot everything from procId == 0
	if (procId == 0) {
		// Get the data provider and give it the points
		_surfacePlot->getDataProvider()->setDataPoints(myPoints);

		// Change the title of the plot and the name of the data
		std::stringstream title;
		// TODO: get the name or comp of the cluster
		title << "First Cluster";
		_surfacePlot->getDataProvider()->setDataName(title.str());
		title << " concentration";
		_surfacePlot->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		_surfacePlot->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
					  << "s";
		_surfacePlot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "surface_TS" << timestep << ".png";
		_surfacePlot->write(fileName.str());
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
} // namespace monitor
} // namespace solver
} // namespace xolotl
