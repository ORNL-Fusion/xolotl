#include <petsc.h>

#include <cassert>

#include <xolotl/factory/interface/MaterialSubOptionsFactory.h>
#include <xolotl/interface/IMaterialSubOptions.h>
#include <xolotl/interface/MultiXolotl.h>
#include <xolotl/interface/XolotlInterface.h>
#include <xolotl/util/GrowthFactorStepSequence.h>
#include <xolotl/util/LinearStepSequence.h>

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

MultiXolotl::MultiXolotl(const std::shared_ptr<ComputeContext>& context,
	const std::shared_ptr<options::IOptions>& options) :
	_computeContext(context),
	_options(options),
	_timeStepper(std::make_unique<util::GrowthFactorStepSequence>(
					 _options->getInitialTimeStep(), _options->getMaxTimeStep(),
					 _options->getTimeStepGrowthFactor()),
		_options->getStartTime(), _options->getEndTime(),
		_options->getNumberOfTimeSteps()),
	_petscContext(std::make_unique<PetscContext>())
{
	// Create primary (whole) network interface
	auto primaryOpts = _options->makeCopy();
	primaryOpts->addProcess("noSolve");
	_primaryInstance = std::make_unique<XolotlInterface>(context, primaryOpts);

	auto subOptions = factory::interface::MaterialSubOptionsFactory::get()
						  .generate(*_options)
						  ->getSubOptions();

	std::vector<std::vector<std::vector<std::uint32_t>>> allBounds;
	std::vector<std::vector<std::vector<xolotl::IdType>>> allMomIdInfo;
	for (auto&& subOpts : subOptions) {
		// Create subinstances
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
		currentTime(), fullConc, std::max((int)temperatures.size() - 2, 1));
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
		currentTime(), fullConc, std::max((int)temperatures.size() - 2, 1));

	// Solve
	for (auto&& sub : _subInstances) {
		// Set the time we want to reach
		sub->setTimes(currentTime(), currentDt());
		// Run the solver
		sub->solveXolotl();
	}
}
} // namespace interface
} // namespace xolotl
