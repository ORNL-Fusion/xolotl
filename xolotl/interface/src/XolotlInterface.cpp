// Includes
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
XolotlInterface::initializeXolotl(int& argc, const char* argv[], MPI_Comm comm)
try {
	computeContext = std::make_unique<ComputeContext>(argc, argv);

	// Initialize the MPI communicator to use
	util::setMPIComm(comm);

	options = std::make_shared<options::Options>();
	options->readParams(argc, argv);

	initializeXolotl();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::initializeXolotl()
try {
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
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::initializeSolver()
try {
	// Initialize the solver
	solver->initialize();
	solverInitialized = true;
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setTimes(double finalTime, double dt)
try {
	// Set the time in the solver
	solver->setTimes(finalTime, dt);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::solveXolotl()
try {
	// Launch the PetscSolver
	solver->solve();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

std::vector<std::vector<std::vector<std::array<double, 4>>>>
XolotlInterface::getLocalNE()
try {
	// Get the solver handler and return the rate vector
	return solverCast(solver)->getSolverHandler()->getLocalNE();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setLocalNE(
	const std::vector<std::vector<std::vector<std::array<double, 4>>>>&
		rateVector)
try {
	// Set the rate vector
	solverCast(solver)->getSolverHandler()->setLocalNE(rateVector);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::getLocalCoordinates(IdType& xs, IdType& xm, IdType& Mx,
	IdType& ys, IdType& ym, IdType& My, IdType& zs, IdType& zm, IdType& Mz)
try {
	// Get the local coordinates
	solverCast(solver)->getSolverHandler()->getLocalCoordinates(
		xs, xm, Mx, ys, ym, My, zs, zm, Mz);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setGBLocation(IdType i, IdType j, IdType k)
try {
	// Set the coordinate of the GB
	solverCast(solver)->getSolverHandler()->setGBLocation(i, j, k);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::resetGBVector()
try {
	// Reset the location
	solverCast(solver)->getSolverHandler()->resetGBVector();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

std::vector<std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
XolotlInterface::getConcVector()
try {
	// Get the vector
	return solver->getConcVector();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setConcVector(std::vector<
	std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
		concVector)
try {
	// Set the vector
	solver->setConcVector(concVector);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

double
XolotlInterface::getPreviousTime()
try {
	return solverCast(solver)->getSolverHandler()->getPreviousTime();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setPreviousTime(double time)
try {
	// Update the fluence from here
	solverCast(solver)->getSolverHandler()->setPreviousTime(time, true);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

double
XolotlInterface::getCurrentDt()
try {
	return solver->getCurrentDt();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setCurrentTimes(double time, double dt)
try {
	solver->setCurrentTimes(time, dt);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

double
XolotlInterface::getNXeGB()
try {
	return solverCast(solver)->getSolverHandler()->getNXeGB();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setNXeGB(double nXe)
try {
	solverCast(solver)->getSolverHandler()->setNXeGB(nXe);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

TS&
XolotlInterface::getTS()
{
	return solver->getTS();
}

std::vector<double>
XolotlInterface::getGridInfo(double& hy, double& hz)
try {
	// Get the solver handler
	auto solverHandler = solverCast(solver)->getSolverHandler();

	// Get the step size
	hy = solverHandler->getStepSizeY();
	hz = solverHandler->getStepSizeZ();

	return solverHandler->getXGrid();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

std::vector<std::vector<AmountType>>
XolotlInterface::getAllClusterBounds()
try {
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();

	return network.getAllClusterBounds();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

std::vector<std::vector<IdType>>
XolotlInterface::getAllMomentIdInfo()
try {
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();

	return network.getAllMomentIdInfo();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::initializeClusterMaps(
	std::vector<std::vector<std::vector<AmountType>>> bounds,
	std::vector<std::vector<std::vector<IdType>>> momIdInfo)
try {
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
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::initializeReactions()
try {
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	network.initializeReactions();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

std::vector<std::vector<std::pair<IdType, double>>>
XolotlInterface::getImplantedFlux()
try {
	// Get the flux handler
	auto fluxHandler = solverCast(solver)->getSolverHandler()->getFluxHandler();

	// Loop on the sub network maps
	std::vector<std::vector<std::pair<IdType, double>>> toReturn;
	for (auto subMap : fromSubNetwork) {
		toReturn.push_back(fluxHandler->getImplantedFlux(subMap));
	}

	return toReturn;
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setImplantedFlux(
	std::vector<std::pair<IdType, double>> fluxVector)
try {
	// Get the flux handler
	auto fluxHandler = solverCast(solver)->getSolverHandler()->getFluxHandler();
	fluxHandler->setImplantedFlux(fluxVector);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setConstantRates(std::vector<std::vector<double>> rates)
try {
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	network.setConstantRates(rates);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

std::vector<std::vector<std::vector<double>>>
XolotlInterface::computeConstantRates(std::vector<std::vector<double>> conc)
try {
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	const auto dof = network.getDOF();

	// Construct the full concentration vector first
	std::vector<double> fullConc(dof, 0.0);
	for (auto i = 0; i < conc.size(); i++)
		for (auto j = 0; j < conc[i].size(); j++) {
			fullConc[fromSubNetwork[i][j]] = conc[i][j];
		}

	auto hConcs =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
			fullConc.data(), dof);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof);
	deep_copy(dConcs, hConcs);

	// Loop on the sub network maps
	std::vector<std::vector<std::vector<double>>> toReturn;
	for (auto l = 0; l < fromSubNetwork.size(); l++) {
		// Get the sub DOF and initialize the rate map
		auto subDOF = fromSubNetwork[l].size();
		std::vector<std::vector<double>> rateMap =
			std::vector(subDOF, std::vector(subDOF + 1, 0.0));
		auto dRates = Kokkos::View<double**>("dRates", subDOF, subDOF + 1);
		auto hRates = Kokkos::create_mirror_view(dRates);
		network.computeConstantRates(dConcs, dRates, l);

		deep_copy(hRates, dRates);
		// Copy element by element
		for (auto i = 0; i < rateMap.size(); i++)
			for (auto j = 0; j < rateMap[0].size(); j++) {
				rateMap[i][j] = hRates(i, j);
			}
		toReturn.push_back(rateMap);
	}

	return toReturn;
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

std::vector<std::vector<std::vector<bool>>>
XolotlInterface::getConstantConnectivities()
try {
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	const auto dof = network.getDOF();

	// Loop on the sub network maps
	std::vector<std::vector<std::vector<bool>>> toReturn;
	for (auto l = 0; l < fromSubNetwork.size(); l++) {
		// Get the sub DOF and initialize the connectivity map
		auto subDOF = fromSubNetwork[l].size();
		std::vector<std::vector<bool>> connMap =
			std::vector(subDOF, std::vector(subDOF + 1, false));
		auto dConns = Kokkos::View<bool**>("dConns", subDOF, subDOF + 1);
		auto hConns = Kokkos::create_mirror_view(dConns);
		network.getConstantConnectivities(dConns, l);

		deep_copy(hConns, dConns);
		// Copy element by element
		for (auto i = 0; i < connMap.size(); i++)
			for (auto j = 0; j < connMap[0].size(); j++) {
				connMap[i][j] = hConns(i, j);
			}
		toReturn.push_back(connMap);
	}

	return toReturn;
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setConstantConnectivities(std::vector<std::vector<bool>> conns)
try {
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	network.setConstantConnectivities(conns);

	return;
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::outputData(double time, std::vector<std::vector<double>> conc)
try {
	// Get the network
	auto& network = solverCast(solver)->getSolverHandler()->getNetwork();
	const auto dof = network.getDOF();
	auto networkSize = network.getNumClusters();

	if (time == 0.0) {
		// Create/open the output files
		std::fstream outputFile;
		outputFile.open("FullAlphaZr.dat", std::fstream::out);
		outputFile << "#time ";
		outputFile << network.getHeaderString();
		outputFile << std::endl;
		outputFile.close();
	}

	// Construct the full concentration vector first
	std::vector<double> fullConc(dof, 0.0);
	for (auto i = 0; i < conc.size(); i++)
		for (auto j = 0; j < conc[i].size(); j++) {
			fullConc[fromSubNetwork[i][j]] = conc[i][j];
		}

	// Set the output precision
	const int outputPrecision = 5;

	// Open the output file
	std::fstream outputFile;
	outputFile.open("FullAlphaZr.dat", std::fstream::out | std::fstream::app);
	outputFile << std::setprecision(outputPrecision);

	auto numSpecies = network.getSpeciesListSize();
	auto myData = std::vector<double>(numSpecies * 6, 0.0);

	// Get the minimum size for the loop densities and diameters
	auto minSizes = solverCast(solver)->getSolverHandler()->getMinSizes();

	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(fullConc.data(), dof);
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

	// Output the data
	outputFile << time << " ";
	for (auto i = 0; i < numSpecies; ++i) {
		outputFile << myData[i * 6] << " " << myData[(i * 6) + 1] << " "
				   << myData[(i * 6) + 2] << " " << myData[(i * 6) + 3] << " "
				   << myData[(i * 6) + 4] << " " << myData[(i * 6) + 5] << " ";
	}
	outputFile << std::endl;

	// Close the output file
	outputFile.close();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

bool
XolotlInterface::getConvergenceStatus()
try {
	return solver->getConvergenceStatus();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::finalizeXolotl()
try {
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
catch (const std::exception& e) {
	reportException(e);
	throw;
}

} /* namespace interface */
} /* namespace xolotl */
