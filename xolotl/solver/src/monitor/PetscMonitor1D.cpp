#include <xolotl/core/Constants.h>
#include <xolotl/core/network/AlloyReactionNetwork.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/io/XFile.h>
#include <xolotl/perf/ScopedTimer.h>
#include <xolotl/solver/PetscSolver.h>
#include <xolotl/solver/monitor/PetscMonitor1D.h>
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
PetscErrorCode
checkNegative(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<PetscMonitor1D*>(ictx)->checkNegative(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
computeTRIDYN(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<PetscMonitor1D*>(ictx)->computeTRIDYN(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
profileTemperature(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr =
		static_cast<PetscMonitor1D*>(ictx)->profileTemperature(
			ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
monitorSeries(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<PetscMonitor1D*>(ictx)->monitorSeries(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

void
PetscMonitor1D::setup(int loop)
{
	PetscErrorCode ierr;

	_loopNumber = loop;

	auto perfHandler = _solverHandler->getPerfHandler();

	// Initialize the timers, including the one for this function.
	_initTimer = perfHandler->getTimer("monitor1D:init");
	perf::ScopedTimer myTimer(_initTimer);
	_checkNegativeTimer = perfHandler->getTimer("monitor1D:checkNeg");
	_tridynTimer = perfHandler->getTimer("monitor1D:tridyn");
	_startStopTimer = perfHandler->getTimer("monitor1D:startStop");
	_heRetentionTimer = perfHandler->getTimer("monitor1D:heRet");
	_xeRetentionTimer = perfHandler->getTimer("monitor1D:xeRet");
	_scatterTimer = perfHandler->getTimer("monitor1D:scatter");
	_seriesTimer = perfHandler->getTimer("monitor1D:series");
	_eventFuncTimer = perfHandler->getTimer("monitor1D:event");
	_postEventFuncTimer = perfHandler->getTimer("monitor1D:postEvent");

	// Get the process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get xolotlViz handler registry
	auto vizHandlerRegistry = _solverHandler->getVizHandler();

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

	// Get the network and its size
	auto& network = _solverHandler->getNetwork();
	const auto networkSize = network.getNumClusters();
	// Get the number of species
	auto numSpecies = network.getSpeciesListSize();

	// Create data depending on the boundary conditions
	if (_solverHandler->getLeftOffset() == 1) {
		_nSurf = std::vector<double>(numSpecies, 0.0);
		_previousSurfFlux = std::vector<double>(numSpecies, 0.0);
	}
	if (_solverHandler->getRightOffset() == 1) {
		_nBulk = std::vector<double>(numSpecies, 0.0);
		_previousBulkFlux = std::vector<double>(numSpecies, 0.0);
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
		if (hasConcentrations and _loopNumber == 0) {
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
		ierr = TSSetPostStep(_ts, checkTimeStep);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSSetPostStep (checkTimeStep) failed.");
	}

	// Set the monitor to check the negative concentrations
	if (flagNeg) {
		// Find the stride to know how often we want to check
		PetscBool flag;
		ierr = PetscOptionsGetReal(
			NULL, NULL, "-check_negative", &_negThreshold, &flag);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsGetReal (-check_negative) "
			"failed.");
		if (!flag)
			_negThreshold = 1.0e-30;

		// checkNegative1D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::checkNegative, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (checkNegative) failed.");
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
		if (hasConcentrations and _loopNumber == 0) {
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
				_nSurf[i] = lastTsGroup->readData1D(nName.str());

				// Create the previous flux attribute name
				std::ostringstream prevFluxName;
				prevFluxName << "previousFlux" << names[i] << "Surf";
				// Read the attribute
				_previousSurfFlux[i] =
					lastTsGroup->readData1D(prevFluxName.str());
			}

			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			_solverHandler->setPreviousTime(previousTime);
		}

		// Get the sputtering yield
		_sputteringYield = _solverHandler->getSputteringYield();

		// Master process
		if (procId == 0 and _loopNumber == 0) {
			// Clear the file where the surface will be written
			std::ofstream outputFile;
			outputFile.open("surface.txt");
			outputFile << "#time height" << std::endl;
			outputFile.close();
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
			"setupPetsc1DMonitor: TSSetEventHandler (eventFunction1D) failed.");
	}

	// Set the monitor to save 1D plot of xenon distribution
	if (flag1DPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a ScatterPlot
			_scatterPlot = vizHandlerRegistry->getPlot(viz::PlotType::SCATTER);

			_scatterPlot->setLogScale();

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
		}

		// monitorScatter1D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorScatter, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (monitorScatter) failed.");
	}

	// Set the monitor to save 1D plot of many concentrations
	if (flagSeries) {
		// Only the master process will create the plot
		if (procId == 0) {
			// Create a ScatterPlot
			_seriesPlot = vizHandlerRegistry->getPlot(viz::PlotType::SERIES);

			// set the log scale
			//			_seriesPlot->setLogScale();

			// Create and set the label provider
			auto labelProvider = std::make_shared<viz::LabelProvider>();
			labelProvider->axis1Label = "x Position on the Grid";
			labelProvider->axis2Label = "Concentration";

			// Give it to the plot
			_seriesPlot->setLabelProvider(labelProvider);

			// To plot a maximum of 18 clusters of the whole benchmark
			const auto loopSize = std::min(18, (int)networkSize);

			// Create a data provider for each cluster in the network
			for (auto i = 0; i < loopSize; i++) {
				// Create the data provider
				auto dataProvider =
					std::make_shared<viz::dataprovider::CvsXDataProvider>();

				// Give it to the plot
				_seriesPlot->addDataProvider(dataProvider);
			}
		}

		// monitorSeries1D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorSeries, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (monitorSeries) failed.");
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
			ierr, "setupPetsc1DMonitor: TSMonitorSet (monitorPerf) failed.");
	}

	// Set the monitor to compute the helium retention
	if (flagHeRetention) {
		auto fluxHandler = _solverHandler->getFluxHandler();
		// Get the previous time if concentrations were stored and initialize
		// the fluence
		if (hasConcentrations and _loopNumber == 0) {
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
					_nSurf[i] = lastTsGroup->readData1D(nName.str());

					// Create the previous flux attribute name
					std::ostringstream prevFluxName;
					prevFluxName << "previousFlux" << names[i] << "Surf";
					// Read the attribute
					_previousSurfFlux[i] =
						lastTsGroup->readData1D(prevFluxName.str());
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
					_nBulk[i] = lastTsGroup->readData1D(nName.str());

					// Create the previous flux attribute name
					std::ostringstream prevFluxName;
					prevFluxName << "previousFlux" << names[i] << "Bulk";
					// Read the attribute
					_previousBulkFlux[i] =
						lastTsGroup->readData1D(prevFluxName.str());
				}
			}
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeFluence, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeHeliumRetention1D will be called at each timestep
		ierr =
			TSMonitorSet(_ts, monitor::computeHeliumRetention, this, nullptr);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (computeHeliumRetention) "
			"failed.");

		// Master process
		if (procId == 0 and _loopNumber == 0) {
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
			if (_solverHandler->largeBubble())
				outputFile << "C_b av_He av_V" << std::endl;
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
		checkPetscError(ierr, "setupPetsc1DMonitor: TSGetDM failed.");
		// Get the local boundaries
		PetscInt xm;
		ierr = DMDAGetCorners(da, NULL, NULL, NULL, &xm, NULL, NULL);
		checkPetscError(ierr, "setupPetsc1DMonitor: DMDAGetCorners failed.");
		// Create the local vectors on each process
		_solverHandler->createLocalNE(xm);

		// Get the previous time if concentrations were stored and initialize
		// the fluence
		if (hasConcentrations and _loopNumber == 0) {
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
			ierr, "setupPetsc1DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeXenonRetention1D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeXenonRetention, this, nullptr);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (computeXenonRetention) "
			"failed.");

		// Master process
		if (procId == 0 and _loopNumber == 0) {
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("retentionOut.txt");
			outputFile << "#time Xenon_content radius partial_radius "
						  "partial_bubble_conc partial_size"
					   << std::endl;
			outputFile.close();
		}
	}

	// Set the monitor to output data for TRIDYN
	if (flagTRIDYN) {
		// computeTRIDYN will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeTRIDYN, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (computeTRIDYN) failed.");
	}

	// Set the monitor to output data for Alloy
	if (flagAlloy) {
		if (procId == 0 and _loopNumber == 0) {
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
		}

		// computeAlloy1D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::computeAlloy, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (computeAlloy) failed.");
	}

	// Set the monitor to compute the temperature profile
	if (flagTemp) {
		if (procId == 0 and _loopNumber == 0) {
			// Uncomment to clear the file where the retention will be written
			std::ofstream outputFile;
			outputFile.open("tempProf.txt");

			// Get the da from _ts
			DM da;
			ierr = TSGetDM(_ts, &da);
			checkPetscError(ierr, "setupPetsc1DMonitor: TSGetDM failed.");

			// Get the total size of the grid
			PetscInt Mx;
			ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE,
				PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
				PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
				PETSC_IGNORE, PETSC_IGNORE);
			checkPetscError(ierr, "setupPetsc1DMonitor: DMDAGetInfo failed.");

			// Get the physical grid
			auto grid = _solverHandler->getXGrid();

			// Loop on the entire grid
			for (auto xi = _solverHandler->getLeftOffset();
				 xi < Mx - _solverHandler->getRightOffset(); xi++) {
				// Set x
				double x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];
				outputFile << x << " ";
			}
			outputFile << std::endl;
			outputFile.close();
		}

		// computeCumulativeHelium1D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::profileTemperature, this, nullptr);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: TSMonitorSet (profileTemperature) failed.");
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
			"setupPetsc1DMonitor: PetscOptionsGetReal (-largest_conc) failed.");

		// monitorLargest1D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::monitorLargest, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (monitorLargest) failed.");
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr =
			PetscOptionsGetReal(NULL, NULL, "-start_stop", &_hdf5Stride, &flag);
		checkPetscError(ierr,
			"setupPetsc1DMonitor: PetscOptionsGetReal (-start_stop) failed.");
		if (!flag)
			_hdf5Stride = 1.0;

		// Compute the correct _hdf5Previous for a restart
		// Get the last time step written in the HDF5 file
		if (hasConcentrations and _loopNumber == 0) {
			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double previousTime = lastTsGroup->readPreviousTime();
			_solverHandler->setPreviousTime(previousTime);
			_hdf5Previous = (PetscInt)(previousTime / _hdf5Stride);
		}

		// Don't do anything if both files have the same name
		if (_hdf5OutputName != _solverHandler->getNetworkName() and
			_loopNumber == 0) {
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
			// copy with HDF5"s H5Ocopy implementation than it is
			// when all processes call the copy function.
			// The checkpoint file must be closed before doing this.
			writeNetwork(
				xolotlComm, _hdf5OutputName, _solverHandler->getNetworkName());
		}

		// startStop1D will be called at each timestep
		ierr = TSMonitorSet(_ts, monitor::startStop, this, nullptr);
		checkPetscError(
			ierr, "setupPetsc1DMonitor: TSMonitorSet (startStop) failed.");
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(_ts, monitor::monitorTime, this, nullptr);
	checkPetscError(
		ierr, "setupPetsc1DMonitor: TSMonitorSet (monitorTime) failed.");
}

PetscErrorCode
PetscMonitor1D::monitorLargest(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
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
	for (auto i = xs; i < xs + xm; i++) {
		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[i];
		// Check the concentration
		if (gridPointSolution[_largestClusterId] > _largestThreshold) {
			ierr = TSSetConvergedReason(ts, TS_CONVERGED_USER);
			CHKERRQ(ierr);
			// Send an error
			throw std::runtime_error(
				"\nxolotlSolver::Monitor1D: The largest cluster "
				"concentration is too high!!");
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor1D::startStop(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declaration
	PetscErrorCode ierr;
	const double **solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	perf::ScopedTimer myTimer(_startStopTimer);

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Compute the dt
	double previousTime = _solverHandler->getPreviousTime();
	double dt = time - previousTime;

	// Don't do anything if it is not on the stride
	if (((PetscInt)((time + dt / 10.0) / _hdf5Stride) <= _hdf5Previous) &&
		timestep > 0) {
		PetscFunctionReturn(0);
	}

	// Update the previous time
	if ((PetscInt)((time + dt / 10.0) / _hdf5Stride) > _hdf5Previous)
		_hdf5Previous++;

	// Gets MPI comm
	auto xolotlComm = util::getMPIComm();

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

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();
	// Write it in the file
	tsGroup->writeGrid(grid);

	// Get the names of the species in the network
	auto numSpecies = network.getSpeciesListSize();
	std::vector<std::string> names;
	for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
		names.push_back(network.getSpeciesName(id));
	}

	if (_solverHandler->moveSurface() || _solverHandler->getLeftOffset() == 1) {
		// Write the surface positions and the associated interstitial
		// quantities in the concentration sub group
		tsGroup->writeSurface1D(_nSurf, _previousSurfFlux, names);
	}

	// Write the bottom impurity information if the bottom is a free surface
	if (_solverHandler->getRightOffset() == 1)
		tsGroup->writeBottom1D(_nBulk, _previousBulkFlux, names);

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

	ierr = computeTRIDYN(ts, timestep, time, solution);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor1D::computeHeliumRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	perf::ScopedTimer myTimer(_heRetentionTimer);

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the flux handler that will be used to know the fluence
	auto fluxHandler = _solverHandler->getFluxHandler();
	// Get the diffusion handler
	auto diffusionHandler = _solverHandler->getDiffusionHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();

	// Get the network
	using NetworkType = core::network::IPSIReactionNetwork;
	using AmountType = NetworkType::AmountType;
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	const auto dof = network.getDOF();

	// Get the array of concentration
	PetscReal** solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Check if we are using the large bubble model
	auto largeBubble = _solverHandler->largeBubble();
	// Store the concentration over the grid
	auto numSpecies = network.getSpeciesListSize();
	auto specIdI = network.getInterstitialSpeciesId();
	auto myConcData = std::vector<double>(numSpecies + 3 * largeBubble, 0.0);

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal* gridPointSolution;

	// Loop on the grid
	for (auto xi = xs; xi < xs + xm; xi++) {
		// Boundary conditions
		if (xi < _solverHandler->getLeftOffset() ||
			xi >= Mx - _solverHandler->getRightOffset())
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
		using Quant = core::network::IReactionNetwork::TotalQuantity;
		std::vector<Quant> quant;
		quant.reserve(numSpecies);
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			quant.push_back({Quant::Type::atom, id, 1});
		}
		auto totals = network.getTotalsVec(dConcs, quant);
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			myConcData[id()] += totals[id()] * hx;
		}
		if (largeBubble) {
			myConcData[numSpecies] +=
				gridPointSolution[dof - 3] * hx; // Bubble conc
			myConcData[numSpecies + 1] +=
				gridPointSolution[dof - 2] * hx; // He conc
			myConcData[numSpecies + 2] +=
				gridPointSolution[dof - 1] * hx; // V conc
		}
	}

	// Get the current process ID
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Determine total concentrations for He, D, T.
	auto totalConcData = std::vector<double>(numSpecies + 3 * largeBubble, 0.0);

	MPI_Reduce(myConcData.data(), totalConcData.data(),
		numSpecies + 3 * largeBubble, MPI_DOUBLE, MPI_SUM, 0, xolotlComm);

	// Get the delta time from the previous timestep to this timestep
	double previousTime = _solverHandler->getPreviousTime();
	double dt = time - previousTime;

	// Look at the fluxes leaving the free surface
	if (_solverHandler->getLeftOffset() == 1) {
		// Set the surface position
		auto xi = 1;

		// Value to know on which processor is the surface
		int surfaceProc = 0;

		// Check we are on the right proc
		if (xi >= xs && xi < xs + xm) {
			// Compute the total number of impurities that left at the surface
			if (timestep > 0) {
				for (auto i = 0; i < numSpecies; ++i) {
					if (i == specIdI() && _solverHandler->moveSurface())
						continue;
					_nSurf[i] += _previousSurfFlux[i] * dt;
				}
			}
			auto myFluxData = std::vector<double>(numSpecies, 0.0);

			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[xi];

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

			// Get the vector of diffusing clusters
			auto diffusingIds = diffusionHandler->getDiffusingIds();

			network.updateOutgoingDiffFluxes(
				gridPointSolution, factor, diffusingIds, myFluxData, xi - xs);

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
				3.0 /
					(core::kBoltzmann * distance * distance * distance *
						distance),
				advecClusters, sinkStrengths, myFluxData, xi - xs);

			for (auto i = 0; i < numSpecies; ++i) {
				if (i == specIdI() && _solverHandler->moveSurface())
					continue;
				_previousSurfFlux[i] = myFluxData[i];
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
		for (auto i = 0; i < numSpecies; ++i) {
			countFluxData.push_back(_nSurf[i]);
			countFluxData.push_back(_previousSurfFlux[i]);
		}
		MPI_Bcast(countFluxData.data(), countFluxData.size(), MPI_DOUBLE,
			surfaceId, xolotlComm);

		// Extract impurity data from broadcast buffer.
		for (auto i = 0; i < numSpecies; ++i) {
			_nSurf[i] = countFluxData[2 * i];
			_previousSurfFlux[i] = countFluxData[(2 * i) + 1];
		}
	}

	// Look at the fluxes going in the bulk if the bottom is a free surface
	if (_solverHandler->getRightOffset() == 1) {
		// Set the bottom surface position
		auto xi = Mx - 2;

		// Value to know on which processor is the bottom
		int bottomProc = 0;

		// Check we are on the right proc
		if (xi >= xs && xi < xs + xm) {
			// Compute the total number of impurities that went in the bulk
			if (timestep > 0) {
				for (auto i = 0; i < numSpecies; ++i) {
					_nBulk[i] += _previousBulkFlux[i] * dt;
				}
			}
			auto myFluxData = std::vector<double>(numSpecies, 0.0);

			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[xi];

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

			// Get the vector of diffusing clusters
			auto diffusingIds = diffusionHandler->getDiffusingIds();

			network.updateOutgoingDiffFluxes(
				gridPointSolution, factor, diffusingIds, myFluxData, xi - xs);

			for (auto i = 0; i < numSpecies; ++i) {
				_previousBulkFlux[i] = myFluxData[i];
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
		for (auto i = 0; i < numSpecies; ++i) {
			countFluxData.push_back(_nBulk[i]);
			countFluxData.push_back(_previousBulkFlux[i]);
		}
		MPI_Bcast(countFluxData.data(), countFluxData.size(), MPI_DOUBLE,
			bottomId, xolotlComm);

		// Extract impurity data from broadcast buffer.
		for (auto i = 0; i < numSpecies; ++i) {
			_nBulk[i] = countFluxData[2 * i];
			_previousBulkFlux[i] = countFluxData[(2 * i) + 1];
		}
	}

	// Master process
	if (procId == 0) {
		// Get the fluence
		double fluence = fluxHandler->getFluence();

		// Print the result
		util::StringStream ss;
		ss << std::endl << "Time: " << time << std::endl;
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			ss << network.getSpeciesName(id)
			   << " content = " << totalConcData[id()] << std::endl;
		}
		ss << "Fluence = " << fluence << std::endl << std::endl;
		XOLOTL_LOG << ss.str();

		// Write the retention and the fluence in a file
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt", std::ios::app);
		outputFile << time << " " << fluence << " ";
		for (auto i = 0; i < numSpecies; ++i) {
			outputFile << totalConcData[i] << " ";
		}
		if (_solverHandler->getRightOffset() == 1) {
			for (auto i = 0; i < numSpecies; ++i) {
				outputFile << _nBulk[i] << " ";
			}
		}
		if (_solverHandler->getLeftOffset() == 1) {
			for (auto i = 0; i < numSpecies; ++i) {
				outputFile << _nSurf[i] << " ";
			}
		}
		if (largeBubble) {
			outputFile << totalConcData[numSpecies] << " "
					   << totalConcData[numSpecies + 1] /
					totalConcData[numSpecies]
					   << " "
					   << totalConcData[numSpecies + 2] /
					totalConcData[numSpecies];
		}
		outputFile << std::endl;
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
PetscMonitor1D::computeXenonRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	perf::ScopedTimer myTimer(_xeRetentionTimer);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();

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
	auto minSizes = _solverHandler->getMinSizes();

	// Get Xe_1
	Composition xeComp = Composition::zero();
	xeComp[Spec::Xe] = 1;
	auto xeCluster = network.findCluster(xeComp, plsm::HostMemSpace{});
	auto xeId = xeCluster.getId();

	// Loop on the grid
	for (auto xi = xs; xi < xs + xm; xi++) {
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
		using TQ = core::network::IReactionNetwork::TotalQuantity;
		using Q = TQ::Type;
		using TQA = util::Array<TQ, 7>;
		auto id =
			core::network::SpeciesId(Spec::Xe, network.getSpeciesListSize());
		auto ms = static_cast<AmountType>(minSizes[id()]);
		auto totals = network.getTotals(dConcs,
			TQA{TQ{Q::total, id, 1}, TQ{Q::atom, id, 1}, TQ{Q::radius, id, 1},
				TQ{Q::total, id, ms}, TQ{Q::atom, id, ms},
				TQ{Q::radius, id, ms}, TQ{Q::volume, id, ms}});

		bubbleConcentration += totals[0] * hx;
		xeConcentration += totals[1] * hx;
		radii += totals[2] * hx;
		partialBubbleConcentration += totals[3] * hx;
		partialSize += totals[4] * hx;
		partialRadii += totals[5] * hx;

		_solverHandler->setVolumeFraction(totals[6], xi - xs);

		_solverHandler->setMonomerConc(
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
		// Check we are on the right proc
		if (xi >= xs && xi < xs + xm) {
			double previousXeFlux = std::get<1>(localNE[xi - xs][0][0]);
			globalXeFlux += previousXeFlux * (grid[xi + 1] - grid[xi]);
			// Set the amount in the vector we keep
			_solverHandler->setLocalXeRate(previousXeFlux * dt, xi - xs);
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
			_solverHandler->setPreviousXeFlux(myRate[0], xi - xs);
		}
	}

	// Master process
	if (procId == 0) {
		// Get the number of xenon that went to the GB
		double nXenon = _solverHandler->getNXeGB();

		// Print the result
		XOLOTL_LOG << std::endl
				   << "Time: " << time << std::endl
				   << "Xenon concentration = " << totalConcData[0] << std::endl
				   << "Xenon GB = " << nXenon << std::endl
				   << std::endl;

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

PetscErrorCode
PetscMonitor1D::computeAlloy(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Get the MPI comm
	auto xolotlComm = util::getMPIComm();

	// Get the process ID
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the physical grid and its length
	auto grid = _solverHandler->getXGrid();
	auto xSize = grid.size();

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

	// Loop on the grid
	for (auto xi = xs; xi < xs + xm; xi++) {
		// Boundary conditions
		if (xi < _solverHandler->getLeftOffset() ||
			xi == Mx - _solverHandler->getRightOffset())
			continue;

		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[xi];

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
				TQA{TQ{Q::total, id, 1}, TQ{Q::radius, id, 1},
					TQ{Q::total, id, ms}, TQ{Q::radius, id, ms}});

			myData[4 * id()] += totals[0];
			myData[(4 * id()) + 1] += 2.0 * totals[1] / myData[4 * id()];
			myData[(4 * id()) + 2] += totals[2];
			myData[(4 * id()) + 3] += 2.0 * totals[3] / myData[(4 * id()) + 2];
		}
	}

	// Sum all the concentrations through MPI reduce
	auto globalData = std::vector<double>(myData.size(), 0.0);
	MPI_Reduce(myData.data(), globalData.data(), myData.size(), MPI_DOUBLE,
		MPI_SUM, 0, xolotlComm);

	// Average the data
	if (procId == 0) {
		for (auto i = 0; i < numSpecies; ++i) {
			globalData[(4 * i) + 1] /= globalData[4 * i];
			globalData[(4 * i) + 3] /= globalData[(4 * i) + 2];
			globalData[4 * i] /= (grid[Mx] - grid[1]);
			globalData[(4 * i) + 2] /= (grid[Mx] - grid[1]);
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
			outputFile << globalData[i * 4] << " " << globalData[(i * 4) + 1]
					   << " " << globalData[(i * 4) + 2] << " "
					   << globalData[(i * 4) + 3] << " ";
		}
		outputFile << std::endl;

		// Close the output file
		outputFile.close();
	}

	// Restore the PETSC solution array
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor1D::monitorScatter(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	perf::ScopedTimer myTimer(_scatterTimer);

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

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the network and its size
	using NetworkType = core::network::NEReactionNetwork;
	using Spec = typename NetworkType::Species;
	using Region = typename NetworkType::Region;
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	auto networkSize = network.getNumClusters();

	// Get the index of the middle of the grid
	auto ix = Mx / 2;

	// If the middle is on this process
	if (ix >= xs && ix < xs + xm) {
		// Create a DataPoint vector to store the data to give to the data
		// provider for the visualization
		auto myPoints =
			std::make_shared<std::vector<viz::dataprovider::DataPoint>>();

		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[ix];

		for (auto i = 0; i < networkSize; i++) {
			// Create a Point with the concentration[i] as the value
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
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
					  << "s";
		_scatterPlot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "Scatter_TS" << timestep << ".png";
		_scatterPlot->render(fileName.str());
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor1D::eventFunction(
	TS ts, PetscReal time, Vec solution, PetscScalar* fvalue)
{
	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	perf::ScopedTimer myTimer(_eventFuncTimer);

	_depthPositions.clear();
	fvalue[0] = 1.0, fvalue[1] = 1.0;

	PetscInt tsNumber = -1;
	ierr = TSGetStepNumber(ts, &tsNumber);
	CHKERRQ(ierr);

	// Skip if it is the same TS as before
	if (tsNumber == _previousTSNumber)
		PetscFunctionReturn(0);

	// Set the previous TS number
	_previousTSNumber = tsNumber;

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

	// Get the position of the surface
	auto xi = _solverHandler->getLeftOffset();

	// Get the network
	using NetworkType = core::network::IPSIReactionNetwork;
	auto& network = dynamic_cast<NetworkType&>(_solverHandler->getNetwork());
	// Get the number of species
	auto numSpecies = network.getSpeciesListSize();
	auto specIdI = network.getInterstitialSpeciesId();

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();

	// Get the flux handler to know the flux amplitude.
	auto fluxHandler = _solverHandler->getFluxHandler();
	double heliumFluxAmplitude = fluxHandler->getFluxAmplitude();

	// Get the delta time from the previous timestep to this timestep
	double dt = time - _solverHandler->getPreviousTime();

	// Work of the moving surface first
	if (_solverHandler->moveSurface()) {
		// Write the initial surface position
		if (procId == 0 and tsNumber == 0) {
			std::ofstream outputFile;
			outputFile.open("surface.txt", std::ios::app);
			outputFile << time << " " << grid[grid.size() - 2] - grid[1]
					   << std::endl;
			outputFile.close();
		}

		// Value to know on which processor is the location of the surface,
		// for MPI usage
		int surfaceProc = 0;

		// if xi is on this process
		if (xi >= xs && xi < xs + xm) {
			// Get the concentrations
			gridPointSolution = solutionArray[xi];

			if (tsNumber > 0) {
				// Compute the total density of intersitials that escaped from
				// the surface since last timestep using the stored flux
				_nSurf[specIdI()] += _previousSurfFlux[specIdI()] * dt;

				// Remove the sputtering yield since last timestep
				_nSurf[specIdI()] -=
					_sputteringYield * heliumFluxAmplitude * dt;
			}

			// Initialize the value for the flux
			auto myFlux = std::vector<double>(numSpecies, 0.0);

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

			network.updateOutgoingDiffFluxes(
				gridPointSolution, factor, _iClusterIds, myFlux, xi - xs);

			// Update the previous flux
			_previousSurfFlux[specIdI()] = myFlux[specIdI()];

			// Set the surface processor
			surfaceProc = procId;
		}

		// Get which processor will send the information
		int surfaceId = 0;
		MPI_Allreduce(
			&surfaceProc, &surfaceId, 1, MPI_INT, MPI_SUM, xolotlComm);

		// Send the information to the other processes
		MPI_Bcast(&_nSurf[specIdI()], 1, MPI_DOUBLE, surfaceId, xolotlComm);
		MPI_Bcast(&_previousSurfFlux[specIdI()], 1, MPI_DOUBLE, surfaceId,
			xolotlComm);

		// Now that all the processes have the same value of nInterstitials,
		// compare it to the threshold to now if we should move the surface

		// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
		double threshold = core::tungstenDensity * (grid[xi] - grid[xi - 1]);
		if (_nSurf[specIdI()] > threshold) {
			// The surface is moving
			fvalue[0] = 0;
		}

		// Update the threshold for erosion (the cell size is not the same)
		threshold = core::tungstenDensity * (grid[xi + 1] - grid[xi]);
		// Moving the surface back
		if (_nSurf[specIdI()] < -threshold * 0.9) {
			// The surface is moving
			fvalue[1] = 0;
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor1D::postEventFunction(TS ts, PetscInt nevents, PetscInt eventList[],
	PetscReal time, Vec solution, PetscBool)
{
	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	perf::ScopedTimer myTimer(_postEventFuncTimer);

	// Check if the surface has moved
	if (nevents == 0) {
		PetscFunctionReturn(0);
	}

	// Check if both events happened
	if (nevents == 2)
		throw std::runtime_error(
			"\nxolotlSolver::Monitor1D: This is not supposed to "
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

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the network
	auto& network = _solverHandler->getNetwork();
	auto dof = network.getDOF();

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();

	// Get the flux handler to know the flux amplitude.
	auto fluxHandler = _solverHandler->getFluxHandler();
	double heliumFluxAmplitude = fluxHandler->getFluxAmplitude();

	// Get the delta time from the previous timestep to this timestep
	double previousTime = _solverHandler->getPreviousTime();
	double dt = time - previousTime;

	// Now takes care of moving surface
	bool moving = false;
	bool movingUp = false;
	for (auto i = 0; i < nevents; i++) {
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
	auto xi = _solverHandler->getLeftOffset();

	// Get the network
	using NetworkType = core::network::IPSIReactionNetwork;
	using AmountType = NetworkType::AmountType;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);

	auto specIdI = psiNetwork->getInterstitialSpeciesId();

	// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
	double threshold = core::tungstenDensity * (grid[xi] - grid[xi - 1]);

	if (movingUp) {
		int nGridPoints = 0;
		// Move the surface up until it is smaller than the next threshold
		while (_nSurf[specIdI()] > threshold &&
			_solverHandler->getLeftOffset() - 2 >= 0) {
			// Move the surface higher
			nGridPoints++;
			// Update the number of interstitials
			_nSurf[specIdI()] -= threshold;

			if (nGridPoints > 0) {
				_solverHandler->setSurfaceOffset(nGridPoints);
				ierr = TSSetConvergedReason(ts, TS_CONVERGED_USER);
				CHKERRQ(ierr);
				PetscFunctionReturn(0);
			}
		}
	}

	// Moving the surface back
	else {
		int nGridPoints = 0;
		// Move it back as long as the number of interstitials in negative
		while (_nSurf[specIdI()] < 0.0) {
			// Compute the threshold to a deeper grid point
			threshold = core::tungstenDensity * (grid[xi + 1] - grid[xi]);
			// Set all the concentrations to 0.0 at xi = 1
			// if xi is on this process
			if (xi >= xs && xi < xs + xm) {
				// Get the concentrations at xi = 1
				gridPointSolution = solutionArray[xi];
				// Loop on DOF
				for (auto i = 0; i < dof; i++) {
					gridPointSolution[i] = 0.0;
				}
			}

			// Move the surface deeper
			nGridPoints--;
			xi++;
			// Update the number of interstitials
			_nSurf[specIdI()] += threshold;
		}

		_solverHandler->setSurfaceOffset(nGridPoints);
		ierr = TSSetConvergedReason(ts, TS_CONVERGED_USER);
		CHKERRQ(ierr);
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor1D::checkNegative(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	perf::ScopedTimer myTimer(_checkNegativeTimer);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the network and dof
	auto& network = _solverHandler->getNetwork();
	const auto nClusters = network.getNumClusters();

	// Loop on the local grid
	for (auto i = xs; i < xs + xm; i++) {
		// Get the pointer to the beginning of the solution data for this grid
		// point
		gridPointSolution = solutionArray[i]; // Loop on the concentrations
		for (auto l = 0; l < nClusters; l++) {
			if (gridPointSolution[l] < _negThreshold &&
				gridPointSolution[l] > 0.0) {
				gridPointSolution[l] = _negThreshold;
			}
			else if (gridPointSolution[l] > -_negThreshold &&
				gridPointSolution[l] < 0.0) {
				gridPointSolution[l] = -_negThreshold;
			}
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOF(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor1D::computeTRIDYN(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	perf::ScopedTimer myTimer(_tridynTimer);

	// Get the MPI communicator
	auto xolotlComm = util::getMPIComm();

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the network
	auto& network = _solverHandler->getNetwork();
	const auto dof = network.getDOF();
	const auto numSpecies = network.getSpeciesListSize();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();

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

	// Save current concentrations as an HDF5 file.
	//
	// First create the file for parallel file access.
	std::ostringstream tdFileStr;
	tdFileStr << "TRIDYN_" << timestep << ".h5";
	io::HDF5File tdFile(tdFileStr.str(),
		io::HDF5File::AccessMode::CreateOrTruncateIfExists, xolotlComm, true);

	// Define a dataset for concentrations.
	// Everyone must create the dataset with the same shape.
	const auto numValsPerGridpoint = 5 + 2;
	const auto firstIdxToWrite = (_solverHandler->getLeftOffset());
	const auto numGridpointsWithConcs = (Mx - firstIdxToWrite);
	io::HDF5File::SimpleDataSpace<2>::Dimensions concsDsetDims = {
		(hsize_t)numGridpointsWithConcs, numValsPerGridpoint};
	io::HDF5File::SimpleDataSpace<2> concsDsetSpace(concsDsetDims);

	const std::string concsDsetName = "concs";
	io::HDF5File::DataSet<double> concsDset(
		tdFile, concsDsetName, concsDsetSpace);

	// Specify the concentrations we will write.
	// We only consider our own grid points.
	const auto myFirstIdxToWrite =
		std::max((IdType)xs, (IdType)firstIdxToWrite);
	auto myEndIdx = (xs + xm); // "end" in the C++ sense; i.e., one-past-last
	auto myNumPointsToWrite =
		(myEndIdx > myFirstIdxToWrite) ? (myEndIdx - myFirstIdxToWrite) : 0;
	io::HDF5File::DataSet<double>::DataType2D<numValsPerGridpoint> myConcs(
		myNumPointsToWrite);

	// Get the interpolated temperature
	std::vector<double> localTemp;
	for (PetscInt i = -1; i <= (PetscInt)xm; ++i) {
		auto gridPointSolution = solutionArray[(PetscInt)xs + i];

		// Get the temperature
		localTemp.push_back(gridPointSolution[dof]);
	}
	auto networkTemp = _solverHandler->interpolateTemperature(localTemp);

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
			auto currIdx = (PetscInt)xi - myFirstIdxToWrite;
			myConcs[currIdx][0] = x;
			// Get the total concentrations at this grid point
			using Quant = core::network::IReactionNetwork::TotalQuantity;
			std::vector<Quant> quant;
			quant.reserve(numSpecies);
			for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
				quant.push_back({Quant::Type::atom, id, 1});
			}
			auto totals = network.getTotalsVec(dConcs, quant);
			for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
				myConcs[currIdx][id() + 1] += totals[id()];
			}
			myConcs[currIdx][numSpecies + 1] = networkTemp[currIdx];
		}
	}

	// Write the concs dataset in parallel.
	// (We write only our part.)
	concsDset.parWrite2D<numValsPerGridpoint>(
		xolotlComm, myFirstIdxToWrite - firstIdxToWrite, myConcs);

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, localSolution, &solutionArray);
	CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da, &localSolution);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor1D::profileTemperature(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	// Gets the process ID (important when it is running in parallel)
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the network and dof
	auto& network = _solverHandler->getNetwork();
	const auto dof = network.getDOF();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();

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

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal* gridPointSolution;

	// Create the output file
	std::ofstream outputFile;
	if (procId == 0) {
		outputFile.open("tempProf.txt", std::ios::app);
		outputFile << time;
	}

	// Create the local vector of temperature wrt temperatureGrid
	std::vector<double> localTemperature;
	// Loop on the local grid including ghosts
	for (auto i = xs; i < xs + xm + 2; i++) {
		// Get the pointer to the beginning of the solution data for this
		// grid point
		gridPointSolution = solutionArray[(PetscInt)i - 1];

		// Get the local temperature
		localTemperature.push_back(gridPointSolution[dof]);
	}

	// Interpolate
	auto updatedTemperature =
		_solverHandler->interpolateTemperature(localTemperature);

	// Loop on the entire grid
	for (auto xi = _solverHandler->getLeftOffset();
		 xi < Mx - _solverHandler->getRightOffset(); xi++) {
		// Set x
		double x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];

		double localTemp = 0.0;
		// Check if this process is in charge of xi
		if (xi >= xs && xi < xs + xm) {
			// Get the local temperature
			localTemp = updatedTemperature[xi - xs + 1];
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
	ierr = DMDAVecRestoreArrayDOFRead(da, localSolution, &solutionArray);
	CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da, &localSolution);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor1D::monitorSeries(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// Initial declarations
	PetscErrorCode ierr;
	const double **solutionArray, *gridPointSolution;
	IdType xs, xm, Mx, ys, ym, My, zs, zm, Mz;

	PetscFunctionBeginUser;

	perf::ScopedTimer myTimer(_seriesTimer);

	double x = 0.0;

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

	// Get local coordinates
	_solverHandler->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Get the network and its size
	auto& network = _solverHandler->getNetwork();
	const auto networkSize = network.getNumClusters();

	// Get the physical grid
	auto grid = _solverHandler->getXGrid();

	// To plot a maximum of 18 clusters of the whole benchmark
	const auto loopSize = std::min(18, (int)networkSize);

	if (procId == 0) {
		// Create a DataPoint vector to store the data to give to the data
		// provider for the visualization
		std::vector<std::vector<viz::dataprovider::DataPoint>> myPoints(
			loopSize);

		// Loop on the grid
		for (auto xi = xs; xi < xs + xm; xi++) {
			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[xi];

			for (auto i = 0; i < loopSize; i++) {
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
		for (auto i = 1; i < worldSize; i++) {
			// Get the size of the local grid of that process
			int localSize = 0;
			MPI_Recv(
				&localSize, 1, MPI_INT, i, 20, xolotlComm, MPI_STATUS_IGNORE);

			// Loop on their grid
			for (auto k = 0; k < localSize; k++) {
				// Get the position
				MPI_Recv(
					&x, 1, MPI_DOUBLE, i, 21, xolotlComm, MPI_STATUS_IGNORE);

				for (auto j = 0; j < loopSize; j++) {
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

		for (auto i = 0; i < loopSize; i++) {
			// Get the data provider and give it the points
			auto thePoints =
				std::make_shared<std::vector<viz::dataprovider::DataPoint>>(
					myPoints[i]);
			_seriesPlot->getDataProvider(i)->setDataPoints(thePoints);
			// TODO: get the name or comp of the cluster
			_seriesPlot->getDataProvider(i)->setDataName(std::to_string(i));
		}

		// Change the title of the plot
		std::stringstream title;
		title << "Concentrations";
		_seriesPlot->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		_seriesPlot->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
					  << "s";
		_seriesPlot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << "log_series_TS" << timestep << ".ppm";
		_seriesPlot->render(fileName.str());
	}

	else {
		// Send the value of the local grid size to the master process
		MPI_Send(&xm, 1, MPI_DOUBLE, 0, 20, xolotlComm);

		// Loop on the grid
		for (auto xi = xs; xi < xs + xm; xi++) {
			// Dump x
			x = (grid[xi] + grid[xi + 1]) / 2.0 - grid[1];

			// Get the pointer to the beginning of the solution data for this
			// grid point
			gridPointSolution = solutionArray[xi];

			// Send the value of the local position to the master process
			MPI_Send(&x, 1, MPI_DOUBLE, 0, 21, xolotlComm);

			for (auto i = 0; i < loopSize; i++) {
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
} // namespace monitor
} // namespace solver
} // namespace xolotl
