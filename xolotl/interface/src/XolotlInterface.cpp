// Includes
#include <cassert>
#include <ctime>
#include <fstream>
#include <iostream>

#include <xolotl/factory/perf/PerfHandlerFactory.h>
#include <xolotl/factory/solver/SolverFactory.h>
#include <xolotl/factory/viz/VizHandlerFactory.h>
#include <xolotl/interface/ComputeContext.h>
#include <xolotl/interface/XolotlInterface.h>
#include <xolotl/options/Options.h>
#include <xolotl/perf/IPerfHandler.h>
#include <xolotl/solver/Solver.h>
#include <xolotl/solver/handler/ISolverHandler.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/version.h>

namespace xolotl
{
namespace interface
{
std::shared_ptr<solver::Solver>
solverCast(const std::shared_ptr<solver::ISolver>& solver) noexcept
{
	auto ret = std::dynamic_pointer_cast<solver::Solver>(solver);
	assert(ret.get() != nullptr);
	return ret;
}

void
reportException(const std::exception& e)
{
	XOLOTL_LOG_ERR << e.what();
	util::Log::flush();
	std::cerr << "Aborting." << std::endl;
}

#ifdef XOLOTL_INTERFACE_REPORT_ERRORS
#define TRY try
#define CATCH \
	catch (const std::exception& e) \
	{ \
		reportException(e); \
		throw; \
	}
#else
#define TRY
#define CATCH
#endif

XolotlInterface::XolotlInterface() = default;

XolotlInterface::XolotlInterface(int& argc, const char* argv[], MPI_Comm comm)
{
	initializeXolotl(argc, argv, comm);
	initializedHere = true;
}

XolotlInterface::XolotlInterface(const std::shared_ptr<ComputeContext>& context,
	const std::shared_ptr<options::IOptions>& opts, MPI_Comm comm) :
	computeContext(context),
	options(opts)
{
	util::setMPIComm(comm);
	initializeXolotl();
	initializedHere = true;
}

XolotlInterface::~XolotlInterface()
{
	if (initializedHere) {
		finalizeXolotl();
	}
}

void
XolotlInterface::printSomething()
{
	std::cout << "I'm in Xolotl !!!" << std::endl;
}

void
XolotlInterface::initializeXolotl(
	int& argc, const char* argv[], MPI_Comm comm) TRY
{
	computeContext = std::make_unique<ComputeContext>(argc, argv);

	// Initialize the MPI communicator to use
	util::setMPIComm(comm);

	options = options::createOptions(argc, argv);
	options->readParams(argc, argv);

	initializeXolotl();
}
CATCH

void
XolotlInterface::initializeXolotl() TRY
{
	if (util::getMPIRank() == 0) {
		// Print the start message
		XOLOTL_LOG << "Starting Xolotl (" << getExactVersionString() << ")\n";
		// TODO! Print copyright message
		// Print date and time
		std::time_t currentTime = std::time(NULL);
		XOLOTL_LOG << std::asctime(std::localtime(&currentTime)) << std::flush;
	}

	// Setup the solver
	solver = factory::solver::SolverFactory::get().generate(*options);
	assert(solver);

	auto processMap = options->getProcesses();
	if (processMap["noSolve"]) {
		// The flux handler still need to be initialized
		auto fluxHandler =
			solverCast(solver)->getSolverHandler()->getFluxHandler();
		auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
		fluxHandler->initializeFluxHandler(network, 0, std::vector<double>());

		// The network needs the temperature for the rates
		// TODO: this assumes constant temperature and 0D case
		auto temperature = std::vector<double>(1, options->getTempParam());
		auto depths = std::vector<double>(1, 0.0);
		network.setTemperatures(temperature, depths);

		return;
	}
	// If constant reactions, initialize later
	if (processMap["constant"]) {
		return;
	}

	// Initialize the solver
	solver->initialize();
	solverInitialized = true;
}
CATCH

void
XolotlInterface::initializeSolver() TRY
{
	// Initialize the solver
	solver->initialize();
	solverInitialized = true;
}
CATCH

void
XolotlInterface::getNetworkTemperature(
	std::vector<double>& temperatures, std::vector<double>& depths) TRY
{
	// Get the temperature and associated depth from the solver handler
	solverCast(solver)->getSolverHandler()->getNetworkTemperature(
		temperatures, depths);

	return;
}
CATCH

void
XolotlInterface::setNetworkTemperature(
	std::vector<double> temperatures, std::vector<double> depths) TRY
{
	solverCast(solver)->getSolverHandler()->getNetwork().setTemperatures(
		temperatures, depths);

	return;
}
CATCH

void
XolotlInterface::setTimes(double finalTime, double dt) TRY
{
	// Set the time in the solver
	solver->setTimes(finalTime, dt);
}
CATCH

void
XolotlInterface::solveXolotl() TRY
{
	// Launch the PetscSolver
	solver->solve();
}
CATCH

std::vector<std::vector<std::vector<std::array<double, 4>>>>
XolotlInterface::getLocalNE() TRY
{
	// Get the solver handler and return the rate vector
	return solverCast(solver)->getSolverHandler()->getLocalNE();
}
CATCH

void
XolotlInterface::setLocalNE(
	const std::vector<std::vector<std::vector<std::array<double, 4>>>>&
		rateVector) TRY
{
	// Set the rate vector
	solverCast(solver)->getSolverHandler()->setLocalNE(rateVector);
}
CATCH

void
XolotlInterface::getLocalCoordinates(IdType& xs, IdType& xm, IdType& Mx,
	IdType& ys, IdType& ym, IdType& My, IdType& zs, IdType& zm, IdType& Mz) TRY
{
	// Get the local coordinates
	solverCast(solver)->getSolverHandler()->getLocalCoordinates(
		xs, xm, Mx, ys, ym, My, zs, zm, Mz);
}
CATCH

void
XolotlInterface::setGBLocation(IdType i, IdType j, IdType k) TRY
{
	// Set the coordinate of the GB
	solverCast(solver)->getSolverHandler()->setGBLocation(i, j, k);
}
CATCH

void
XolotlInterface::resetGBVector() TRY
{
	// Reset the location
	solverCast(solver)->getSolverHandler()->resetGBVector();
}
CATCH

std::vector<std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
XolotlInterface::getConcVector() TRY
{
	// Get the vector
	return solver->getConcVector();
}
CATCH

void
XolotlInterface::setConcVector(std::vector<
	std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
		concVector) TRY
{
	// Set the vector
	solver->setConcVector(concVector);
}
CATCH

double
XolotlInterface::getPreviousTime() TRY
{
	return solverCast(solver)->getSolverHandler()->getPreviousTime();
}
CATCH

void
XolotlInterface::setPreviousTime(double time) TRY
{
	// Update the fluence from here
	solverCast(solver)->getSolverHandler()->setPreviousTime(time, true);
}
CATCH

double
XolotlInterface::getCurrentDt() TRY
{
	return solver->getCurrentDt();
}
CATCH

void
XolotlInterface::setCurrentTimes(double time, double dt) TRY
{
	solver->setCurrentTimes(time, dt);
}
CATCH

double
XolotlInterface::getNXeGB() TRY
{
	return solverCast(solver)->getSolverHandler()->getNXeGB();
}
CATCH

void
XolotlInterface::setNXeGB(double nXe) TRY
{
	solverCast(solver)->getSolverHandler()->setNXeGB(nXe);
}
CATCH

TS&
XolotlInterface::getTS() TRY
{
	return solver->getTS();
}
CATCH

std::vector<double>
XolotlInterface::getGridInfo(double& hy, double& hz) TRY
{
	// Get the solver handler
	auto solverHandler = solverCast(solver)->getSolverHandler();

	// Get the step size
	hy = solverHandler->getStepSizeY();
	hz = solverHandler->getStepSizeZ();

	return solverHandler->getXGrid();
}
CATCH

std::vector<std::vector<AmountType>>
XolotlInterface::getAllClusterBounds() TRY
{
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();

	return network.getAllClusterBounds();
}
CATCH

std::vector<std::vector<IdType>>
XolotlInterface::getAllMomentIdInfo() TRY
{
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();

	return network.getAllMomentIdInfo();
}
CATCH

void
XolotlInterface::initializeClusterMaps(
	std::vector<std::vector<std::vector<AmountType>>> bounds,
	std::vector<std::vector<std::vector<IdType>>> momIdInfo) TRY
{
	// Create the local maps
	auto currentBounds = getAllClusterBounds();
	auto currentMomIdInfo = getAllMomentIdInfo();
	// Loop on the sub network bounds
	for (auto i = 0; i < bounds.size(); i++) {
		auto subBounds = bounds[i];
		// Create a new vector
		std::vector<IdType> temp;
		// Loop on the cluster entries
		for (auto bound : subBounds) {
			// Look for the same cluster in currentBounds
			for (auto j = 0; j < currentBounds.size(); j++) {
				if (bound == currentBounds[j]) {
					temp.push_back(j);
					break;
				}
			}
		}
		// Loop on the momId entries
		auto subMomIdInfo = momIdInfo[i];
		for (auto l = 0; l < subMomIdInfo.size(); l++) {
			auto idMap = subMomIdInfo[l];
			for (auto j = 0; j < idMap.size(); j++) {
				temp.push_back(currentMomIdInfo[temp[l]][j]);
			}
		}
		fromSubNetwork.push_back(temp);
	}

	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	network.initializeClusterMap(bounds, momIdInfo, fromSubNetwork);
}
CATCH

void
XolotlInterface::initializeReactions() TRY
{
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	network.initializeReactions();
}
CATCH

std::vector<std::vector<std::pair<IdType, double>>>
XolotlInterface::getImplantedFlux() TRY
{
	// Get the flux handler
	auto fluxHandler = solverCast(solver)->getSolverHandler()->getFluxHandler();

	// Loop on the sub network maps
	std::vector<std::vector<std::pair<IdType, double>>> toReturn;
	for (auto subMap : fromSubNetwork) {
		toReturn.push_back(fluxHandler->getImplantedFlux(subMap));
	}

	return toReturn;
}
CATCH

void
XolotlInterface::setImplantedFlux(
	std::vector<std::pair<IdType, double>> fluxVector) TRY
{
	// Get the flux handler
	auto fluxHandler = solverCast(solver)->getSolverHandler()->getFluxHandler();
	fluxHandler->setImplantedFlux(fluxVector);
}
CATCH

struct RatesCapsule
{
	core::network::IReactionNetwork::RatesView view;
};

std::shared_ptr<RatesCapsule>
XolotlInterface::makeRatesCapsule() const TRY
{
	return std::make_shared<RatesCapsule>();
}
CATCH

void
XolotlInterface::setConstantRates(
	const std::shared_ptr<RatesCapsule>& rates, IdType gridIndex) TRY
{
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	network.setConstantRates(rates->view, gridIndex);
}
CATCH

void
XolotlInterface::computeConstantRates(std::vector<std::vector<double>> conc,
	IdType gridIndex, std::vector<std::shared_ptr<RatesCapsule>>& rates) TRY
{
	assert(rates.size() == fromSubNetwork.size());

	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	const auto dof = network.getDOF();

	// Construct the full concentration vector first
	std::vector<double> fullConc(dof, 0.0);
	for (auto i = 0; i < conc.size(); i++) {
		for (auto j = 0; j < conc[i].size(); j++) {
			fullConc[fromSubNetwork[i][j]] = conc[i][j];
		}
	}

	auto hConcs =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
			fullConc.data(), dof);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof);
	deep_copy(dConcs, hConcs);

	// Loop on the sub network maps
	for (auto l = 0; l < fromSubNetwork.size(); l++) {
		// TODO: should this allocation happen only once (in makeRatesCapsule)?
		//       (we'd still need to zero out the data here)
		rates[l]->view = Kokkos::View<double*>("dRates", _subEntries[l]);
		network.computeConstantRates(dConcs, rates[l]->view, l, gridIndex);
	}
}
CATCH

std::vector<std::pair<std::vector<IdType>, std::vector<IdType>>>
XolotlInterface::getConstantConnectivities() TRY
{
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();

	// Loop on the sub network maps
	std::vector<std::pair<std::vector<IdType>, std::vector<IdType>>> toReturn;
	for (auto l = 0; l < fromSubNetwork.size(); l++) {
		// Get the sub DOF and initialize the connectivity map
		auto subDOF = fromSubNetwork[l].size();
		std::vector<IdType> rows;
		std::vector<IdType> entries;
		auto dConns = Kokkos::View<bool**>("dConns", subDOF, subDOF + 1);
		auto hConns = Kokkos::create_mirror_view(dConns);
		network.getConstantConnectivities(dConns, l);

		deep_copy(hConns, dConns);
		// Create the sparse connectivity
		rows.push_back(0);
		for (auto i = 0; i < hConns.extent(0); i++) {
			IdType count = 0;
			for (auto j = 0; j < hConns.extent(1); j++) {
				if (hConns(i, j)) {
					entries.push_back(j);
					count++;
				}
			}
			rows.push_back(rows.back() + count);
		}
		_subEntries.push_back(entries.size());
		toReturn.push_back(
			std::make_pair<std::vector<IdType>, std::vector<IdType>>(
				std::move(rows), std::move(entries)));
	}

	return toReturn;
}
CATCH

void
XolotlInterface::initializeRateEntries(
	const std::vector<std::pair<std::vector<IdType>, std::vector<IdType>>>&
		conns) TRY
{
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	network.initializeRateEntries(conns);
}
CATCH

void
XolotlInterface::setConstantConnectivities(
	std::pair<std::vector<IdType>, std::vector<IdType>> conns) TRY
{
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	network.setConstantConnectivities(conns);

	return;
}
CATCH

void
XolotlInterface::outputData(double time,
	std::vector<std::vector<std::vector<double>>> conc, IdType localSize) TRY
{
	// Get the MPI comm
	auto xolotlComm = util::getMPIComm();

	// Get the process ID
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);

	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	const auto dof = network.getDOF();

	if (time == 0.0 and procId == 0) {
		network.writeMonitorOutputHeader();
	}

	auto numSpecies = network.getSpeciesListSize();
	auto myData = std::vector<double>(network.getMonitorDataLineSize(), 0.0);

	for (auto xi = 0; xi < localSize; xi++) {
		// Construct the full concentration vector first
		std::vector<double> fullConc(dof, 0.0);
		for (auto i = 0; i < conc[xi].size(); i++)
			for (auto j = 0; j < conc[xi][i].size(); j++) {
				fullConc[fromSubNetwork[i][j]] = conc[xi][i][j];
			}

		using HostUnmanaged =
			Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
		auto hConcs = HostUnmanaged(fullConc.data(), dof);
		auto dConcs = Kokkos::View<double*>("Concentrations", dof);
		deep_copy(dConcs, hConcs);

		network.addMonitorDataValues(dConcs, 1.0, myData);
	}

	network.writeMonitorDataLine(myData, time);
}
CATCH

bool
XolotlInterface::getConvergenceStatus() TRY
{
	return solver->getConvergenceStatus();
}
CATCH

void
XolotlInterface::finalizeXolotl() TRY
{
	// Call solver finalize
	if (solverInitialized) {
		solver->finalize();
	}

	auto perfHandler = solverCast(solver)->getSolverHandler()->getPerfHandler();

	// Get the MPI rank
	int rank = util::getMPIRank();

	if (options->usePerfOutputYAML()) {
		auto ofs = std::ofstream("perf_r" + std::to_string(rank) + ".yaml");
		perfHandler->reportData(ofs);
	}

	// Report statistics about the performance data collected during
	// the run we just completed.
	perf::PerfObjStatsMap<perf::ITimer::ValType> timerStats;
	perf::PerfObjStatsMap<perf::IEventCounter::ValType> counterStats;
	perf::PerfObjStatsMap<perf::IHardwareCounter::CounterType> hwCtrStats;
	perfHandler->collectStatistics(timerStats, counterStats, hwCtrStats);

	if (rank == 0) {
		util::StringStream ss;
		perfHandler->reportStatistics(ss, timerStats, counterStats, hwCtrStats);
		XOLOTL_LOG << ss.str();
	}

	solver.reset();
}
CATCH

} /* namespace interface */
} /* namespace xolotl */
