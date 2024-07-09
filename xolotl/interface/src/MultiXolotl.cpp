#include <petsc.h>

#include <cassert>
#include <fstream>

#include <xolotl/factory/interface/MaterialSubOptionsFactory.h>
#include <xolotl/interface/IMaterialSubOptions.h>
#include <xolotl/interface/MultiXolotl.h>
#include <xolotl/interface/XolotlInterface.h>
#include <xolotl/util/GrowthFactorStepSequence.h>
#include <xolotl/util/LinearStepSequence.h>
#include <xolotl/util/Log.h>

namespace xolotl
{
namespace interface
{
struct PetscContext
{
	PetscContext()
	{
		PetscInitialize(nullptr, nullptr, nullptr, nullptr);
	}

	~PetscContext()
	{
		PetscFinalize();
	}
};

auto
readRestartFile(const std::string& fileName)
{
	auto ret = std::vector<std::string>{};
	if (fileName.empty()) {
		return ret;
	}

	auto ifs = std::ifstream(fileName);
	if (!ifs) {
		throw std::runtime_error("Problem opening file to read: " + fileName);
	}

	std::string line;
	while (std::getline(ifs, line)) {
		ret.push_back(line);
	}

	return ret;
}

void
writeRestartFile(const std::vector<std::string>& fileNames)
{
	std::string name = "xolotlRestart.txt";
	auto ofs = std::ofstream(name);
	if (!ofs) {
		throw std::runtime_error("Problem opening file to write: " + name);
	}
	for (auto&& name : fileNames) {
		ofs << name << '\n';
	}
	ofs.close();
	XOLOTL_LOG << "MultiXolotl: Restart file written to: " << name;
}

auto
readStopData()
{
	std::string name = "xolotlStop_data.txt";
	auto ifs = std::ifstream(name);
	if (!ifs) {
		throw std::runtime_error("Problem opening file to read: " + name);
	}
	double time, dt;
	std::size_t step;
	ifs >> step >> time >> dt;
	return std::make_tuple(step, time, dt);
}

void
MultiXolotl::writeStopData()
{
	std::string name = "xolotlStop_data.txt";
	auto ofs = std::ofstream(name);
	if (!ofs) {
		throw std::runtime_error("Problem opening file to write: " + name);
	}
	ofs << currentStep() << ' ' << currentTime() << ' ' << currentDt() << '\n';
}

util::TimeStepper
makeTimeStepper(const options::IOptions& options)
{
	auto initDt = options.getInitialTimeStep();
	auto maxDt = options.getMaxTimeStep();
	auto gf = options.getTimeStepGrowthFactor();
	auto startTime = options.getStartTime();
	auto endTime = options.getEndTime();
	auto maxSteps = options.getNumberOfTimeSteps();

	std::size_t step = 0;

	if (!options.getRestartFilePath().empty()) {
		std::tie(step, startTime, initDt) = readStopData();
	}

	auto sequence = std::make_unique<util::GrowthFactorStepSequence>(
		initDt, maxDt, gf, step);

	return util::TimeStepper(std::move(sequence), startTime, endTime, maxSteps);
}

MultiXolotl::MultiXolotl(const std::shared_ptr<ComputeContext>& context,
	const std::shared_ptr<options::IOptions>& options) :
	_computeContext(context),
	_options(options),
	_timeStepper(makeTimeStepper(*options)),
	_petscContext(std::make_unique<PetscContext>())
{
	// Check for restart
	auto restartFile = _options->getRestartFilePath();
	bool restarting = !restartFile.empty();

	// Check for checkpoint
	_checkpointing =
		(_options->getPetscArg().find("-start_stop") != std::string::npos);
	std::vector<std::string> checkpointFiles{"xolotlStop_0.h5"};

	// Create primary (whole) network interface
	auto primaryOpts = _options->makeCopy();
	primaryOpts->setCheckpointFilePath(checkpointFiles[0]); // don't need one?
	primaryOpts->addProcess("noSolve");
	_primaryInstance = std::make_unique<XolotlInterface>(context, primaryOpts);

	// Generate suboptions
	auto subOptions = factory::interface::MaterialSubOptionsFactory::get()
						  .generate(*_options)
						  ->getSubOptions();

	// Get restart files
	auto restartFiles = readRestartFile(restartFile);
	if (restarting && restartFiles.size() != (subOptions.size() + 1)) {
		auto msgstrm = std::stringstream()
			<< "MultiXolotl: Number of restart files (" << restartFiles.size()
			<< ") must match number of instances (" << (subOptions.size() + 1)
			<< ")";
		XOLOTL_ERROR(std::runtime_error, msgstrm.str());
		// auto msg = msgstrm.str();
		// XOLOTL_LOG_ERR << msg;
		// throw std::runtime_error(msg);
	}
	if (restarting) {
		primaryOpts->setRestartFilePath(restartFiles[0]);
	}

	// Create subinstances
	std::vector<std::vector<std::vector<std::uint32_t>>> allBounds;
	std::vector<std::vector<std::vector<xolotl::IdType>>> allMomIdInfo;
	for (int id = 0; id < subOptions.size(); ++id) {
		auto& subOpts = subOptions[id];

		// Handle restart and checkpoint files
		if (restarting) {
			subOpts->setRestartFilePath(restartFiles[id + 1]);
		}
		const auto& ckFile = checkpointFiles.emplace_back(
			"xolotlStop_" + std::to_string(id + 1) + ".h5");
		subOpts->setCheckpointFilePath(ckFile);

		// Construct subinstance
		subOpts->addProcess("constant");
		auto& sub = _subInstances.emplace_back(
			std::make_unique<XolotlInterface>(context, subOpts));

		// Bounds
		auto bounds = sub->getAllClusterBounds();
		allBounds.push_back(bounds);
		_subDOFs.push_back(bounds.size());

		// Moment Ids
		auto momIdInfo = sub->getAllMomentIdInfo();
		allMomIdInfo.push_back(momIdInfo);

		// Constant rate capsules
		_constantRates.push_back(sub->makeRatesCapsule());
	}

	// Write restart file
	if (_checkpointing) {
		writeRestartFile(checkpointFiles);
	}

	// Pass data to primary
	_primaryInstance->initializeClusterMaps(allBounds, allMomIdInfo);

	// Subnetwork reactions
	auto connectivities = _primaryInstance->getConstantConnectivities();
	assert(connectivities.size() == _subInstances.size());
	_primaryInstance->initializeRateEntries(connectivities);
	for (IdType i = 0; i < _subInstances.size(); ++i) {
		auto& sub = _subInstances[i];
		sub->setConstantConnectivities(connectivities[i]);
		sub->initializeReactions();
		sub->initializeSolver();
	}

	// Fluxes
	auto fluxVector = _primaryInstance->getImplantedFlux();
	for (IdType i = 0; i < _subInstances.size(); ++i) {
		_subInstances[i]->setImplantedFlux(fluxVector[i]);
	}
}

MultiXolotl::~MultiXolotl()
{
	std::vector<double> temperatures;
	std::vector<double> depths;
	_subInstances[0]->getNetworkTemperature(temperatures, depths);
	// Loop on the grid points
	std::vector<std::vector<std::vector<double>>> fullConc;
	// 0D
	if (temperatures.size() < 2) {
		std::vector<std::vector<double>> conc;
		// Loop on the sub interfaces to get all the concentrations
		for (auto i = 0; i < _subInstances.size(); i++) {
			auto sparseConc = _subInstances[i]->getConcVector();
			std::vector<double> subConc(_subDOFs[i], 0.0);
			for (auto pair : sparseConc[0][0][0]) {
				if (pair.first < _subDOFs[i])
					subConc[pair.first] = pair.second;
			}
			conc.push_back(subConc);
		}
		fullConc.push_back(conc);
	}
	// 1D
	else {
		// Loop on the grid points
		for (auto j = 0; j < temperatures.size() - 2; j++) {
			std::vector<std::vector<double>> conc;
			// Loop on the sub interfaces to get all the concentrations
			for (auto i = 0; i < _subInstances.size(); i++) {
				auto sparseConc = _subInstances[i]->getConcVector();
				std::vector<double> subConc(_subDOFs[i], 0.0);
				for (auto pair : sparseConc[0][0][j]) {
					if (pair.first < _subDOFs[i]) {
						subConc[pair.first] = pair.second;
					}
				}
				conc.push_back(subConc);
			}

			fullConc.push_back(conc);
		}
	}

	// Print the result
	_primaryInstance->outputData(
		previousTime(), fullConc, std::max((int)temperatures.size() - 2, 1));

	// Write stop data
	if (_checkpointing) {
		writeStopData();
	}
}

double
MultiXolotl::previousTime() const noexcept
{
	return _timeStepper.previousTime();
}

double
MultiXolotl::currentTime() const noexcept
{
	return _timeStepper.currentTime();
}

double
MultiXolotl::currentDt() const noexcept
{
	return _timeStepper.currentTimeStepSize();
}

std::size_t
MultiXolotl::currentStep() const noexcept
{
	return _timeStepper.currentStep();
}

void
MultiXolotl::solveXolotl()
{
	for (_timeStepper.start(); _timeStepper; ++_timeStepper) {
		solveStep();
	}
}

void
MultiXolotl::solveStep()
{
	// Transfer the temperature to the full network
	std::vector<double> temperatures;
	std::vector<double> depths;
	_subInstances[0]->getNetworkTemperature(temperatures, depths);
	std::vector<std::vector<std::vector<double>>> fullConc;

	// 0D
	if (temperatures.size() < 2) {
		std::vector<std::vector<double>> conc;
		// Loop on the sub interfaces to get all the concentrations
		for (auto i = 0; i < _subInstances.size(); i++) {
			auto sparseConc = _subInstances[i]->getConcVector();
			std::vector<double> subConc(_subDOFs[i], 0.0);
			for (auto pair : sparseConc[0][0][0]) {
				if (pair.first < _subDOFs[i]) {
					subConc[pair.first] = pair.second;
				}
			}
			conc.push_back(subConc);
		}

		fullConc.push_back(conc);

		// Compute the new rates
		std::vector<double> temperature = {temperatures[0]};
		std::vector<double> depth = {depths[0]};
		_primaryInstance->setNetworkTemperature(temperature, depth);
		_primaryInstance->computeConstantRates(conc, 0, _constantRates);

		// Pass them
		for (auto i = 0; i < _subInstances.size(); i++) {
			_subInstances[i]->setConstantRates(_constantRates[i], 0);
		}
	}
	// 1D
	else {
		// Loop on the grid points
		for (auto j = 0; j < temperatures.size() - 2; j++) {
			std::vector<std::vector<double>> conc;
			// Loop on the sub interfaces to get all the concentrations
			for (auto i = 0; i < _subInstances.size(); i++) {
				auto sparseConc = _subInstances[i]->getConcVector();
				std::vector<double> subConc(_subDOFs[i], 0.0);
				for (auto pair : sparseConc[0][0][j]) {
					if (pair.first < _subDOFs[i]) {
						subConc[pair.first] = pair.second;
					}
				}
				conc.push_back(subConc);
			}

			fullConc.push_back(conc);

			// Compute the new rates
			std::vector<double> temperature = {temperatures[j + 1]};
			std::vector<double> depth = {depths[j + 1]};
			_primaryInstance->setNetworkTemperature(temperature, depth);
			_primaryInstance->computeConstantRates(conc, 0, _constantRates);

			// Pass them
			for (auto i = 0; i < _subInstances.size(); i++) {
				_subInstances[i]->setConstantRates(_constantRates[i], j + 1);
			}
		}
	}

	// Print the result
	_primaryInstance->outputData(
		previousTime(), fullConc, std::max((int)temperatures.size() - 2, 1));

	// Solve
	for (auto&& sub : _subInstances) {
		// Set the time we want to reach
		sub->setTimes(currentTime(), currentDt());
		// Provide our current step as the external control step
		sub->setExternalControlStep(currentStep());
		// Run the solver
		sub->solveXolotl();
	}
}
} // namespace interface
} // namespace xolotl
