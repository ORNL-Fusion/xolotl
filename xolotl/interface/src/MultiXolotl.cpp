#include <petsc.h>

#include <xolotl/interface/MultiXolotl.h>
#include <xolotl/interface/XolotlInterface.h>
#include <xolotl/util/GrowthFactorStepSequence.h>
#include <xolotl/util/LinearStepSequence.h>

namespace xolotl
{
namespace interface
{
class PetscContext
{
public:
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
	_petscContext(std::make_unique<PetscContext>()),
	_maxDt(10.0 /* FIXME */)
{
	// Create primary (whole) network interface
	auto primaryOpts = _options->makeCopy();
	primaryOpts->addProcess("noSolve");
	_primaryInstance = std::make_unique<XolotlInterface>(context, primaryOpts);

	std::vector<std::vector<std::vector<std::uint32_t>>> allBounds;
	std::vector<std::vector<std::vector<xolotl::IdType>>> allMomIdInfo;
	const auto& netParams = _options->getNetworkParameters();
	auto tmpParams = netParams;
	for (std::size_t i = 0; i < netParams.size(); ++i) {
		if (netParams[i] == 0) {
			continue;
		}

		// Create subinstances
		tmpParams.assign(netParams.size(), 0);
		tmpParams[i] = netParams[i];
		auto subOpts = _options->makeCopy();
		subOpts->setNetworkParameters(tmpParams);
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
	}

	// Pass data to primary
	_primaryInstance->initializeClusterMaps(allBounds, allMomIdInfo);

	// Subnetwork reactions
	auto connectivities = _primaryInstance->getConstantConnectivities();
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
		for (auto j = 0; j < temperatures.size() - 2; j++) {
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
	}

	// Print the result
	_primaryInstance->outputData(
		_currentTime, fullConc, std::max((int)temperatures.size() - 2, 1));
}

void
MultiXolotl::solveXolotl()
{
	auto seq = util::GrowthFactorStepSequence(1.0, _maxDt, 1.3, 20);

	/*
	 * couplingTimeStepParams: initial, max, growthfactor, maxSteps
	 */

	for (seq.start(); seq; seq.step()) {
		_currentTime = seq.current();
		_currentDt = seq.stepSize();
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
		auto constantRates = _primaryInstance->computeConstantRates(conc, 0);

		// Pass them
		for (auto i = 0; i < _subInstances.size(); i++) {
			_subInstances[i]->setConstantRates(constantRates[i], 0);
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
			auto constantRates =
				_primaryInstance->computeConstantRates(conc, 0);

			// Pass them
			for (auto i = 0; i < _subInstances.size(); i++) {
				_subInstances[i]->setConstantRates(constantRates[i], j + 1);
			}
		}
	}

	// Print the result
	_primaryInstance->outputData(
		_currentTime, fullConc, std::max((int)temperatures.size() - 2, 1));

	// Solve
	for (auto&& sub : _subInstances) {
		// Set the time we want to reach
		sub->setTimes(_currentTime, _currentDt);
		// Run the solver
		sub->solveXolotl();
	}
}

// InputParameters
// XolotlNetworkProblem::validParams()
// {
// 	InputParameters params = ExternalProblem::validParams();

// 	// Parameter for the Xolotl file name
// 	params.addRequiredParam<FileName>("network_xolotl_filename",
// 		"Name with the path for the Xolotl input file with the full network");
// 	params.addRequiredParam<std::vector<FileName>>(
// 		"subnetwork_xolotl_filenames",
// 		"Name with the path for the Xolotl input files");
// 	params.addParam<Real>("max_dt", 1.0e9, "The maximum coupling dt (s)");
// 	return params;
// }

// bool
// XolotlNetworkProblem::converged()
// {
// 	return true;
// }

// void
// XolotlNetworkProblem::saveState()
// {
// 	// Update the values from Xolotl
// 	_conc_vector = _networkInterface->getConcVector();
// 	_current_dt = _networkInterface->getCurrentDt();
// 	_previous_time = _networkInterface->getPreviousTime();

// 	xolotl::IdType i, j, k, xs, ys, zs, xm, ym, zm, Mx, My, Mz;
// 	_networkInterface->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);
// }

// void
// XolotlNetworkProblem::setState()
// {
// 	// Set them in Xolotl
// 	_networkInterface->setConcVector(_conc_vector);
// 	_networkInterface->setCurrentTimes(_current_time, _current_dt);
// 	_networkInterface->setPreviousTime(_previous_time);
// }

// XolotlNetworkProblem::XolotlNetworkProblem(const InputParameters& params) :
// 	ExternalProblem(params),
// 	_network_xolotl_filename(getParam<FileName>("network_xolotl_filename")),
// 	_subnetwork_xolotl_filenames(
// 		getParam<std::vector<FileName>>("subnetwork_xolotl_filenames")),
// 	_current_time(declareRestartableData<Real>("current_time", 0.0)),
// 	_max_dt(getParam<Real>("max_dt")),
// 	_localXM(declareRestartableData<xolotl::IdType>("localXM", 1)),
// 	_current_dt(declareRestartableData<Real>("current_dt", 0.0)),
// 	_previous_time(declareRestartableData<Real>("previous_time", 0.0)),
// 	_conc_vector(declareRestartableData<std::vector<std::vector<
// 			std::vector<std::vector<std::pair<xolotl::IdType, Real>>>>>>(
// 		"conc_vector"))
// {
// 	// Create the whole network interface
// 	int argc = 2;
// 	const char* argv[argc + 1];
// 	std::string fakeAppName = "mainXolotl";
// 	argv[0] = fakeAppName.c_str();
// 	argv[1] = _network_xolotl_filename.c_str();

// 	_networkInterface = std::make_shared<XolotlInterface>();
// 	_networkInterface->initializeXolotl(
// 		argc, argv, (_app.getCommunicator())->get());

// 	_subInterfaces.clear();

// 	// Loop on the number of parameter files
// 	std::vector<std::vector<std::vector<std::uint32_t>>> allBounds;
// 	std::vector<std::vector<std::vector<xolotl::IdType>>> allMomIdInfo;
// 	for (auto name : _subnetwork_xolotl_filenames) {
// 		_subInterfaces.push_back(std::make_shared<XolotlInterface>());

// 		std::string fakeAppName = "subXolotl";
// 		argv[0] = fakeAppName.c_str();
// 		argv[1] = name.c_str();

// 		(_subInterfaces.back())
// 			->initializeXolotl(argc, argv, (_app.getCommunicator())->get());

// 		// Get the bounds
// 		auto bounds = _subInterfaces.back()->getAllClusterBounds();

// 		// Add them to the main vector
// 		allBounds.push_back(bounds);
// 		_subDOFs.push_back(bounds.size());

// 		// Get the mom Id info
// 		auto momIdInfo = _subInterfaces.back()->getAllMomentIdInfo();
// 		allMomIdInfo.push_back(momIdInfo);
// 	}

// 	// Pass it the the network instance
// 	_networkInterface->initializeClusterMaps(allBounds, allMomIdInfo);

// 	// Take care of the sub network reactions
// 	auto connectivities = _networkInterface->getConstantConnectivities();
// 	// Loop on the sub interfaces
// 	xolotl::IdType count = 0;
// 	for (auto inter : _subInterfaces) {
// 		inter->setConstantConnectivities(connectivities[count]);
// 		inter->initializeReactions();
// 		inter->initializeSolver();

// 		count++;
// 	}

// 	// Take care of the fluxes
// 	auto fluxVector = _networkInterface->getImplantedFlux();
// 	for (auto i = 0; i < _subInterfaces.size(); i++) {
// 		_subInterfaces[i]->setImplantedFlux(fluxVector[i]);
// 	}
// }

// XolotlNetworkProblem::~XolotlNetworkProblem()
// {
// 	std::vector<double> temperatures;
// 	std::vector<double> depths;
// 	_subInterfaces[0]->getNetworkTemperature(temperatures, depths);
// 	// Loop on the grid points
// 	std::vector<std::vector<std::vector<double>>> fullConc;
// 	// 0D
// 	if (temperatures.size() < 2) {
// 		std::vector<std::vector<double>> conc;
// 		// Loop on the sub interfaces to get all the concentrations
// 		for (auto i = 0; i < _subInterfaces.size(); i++) {
// 			auto sparseConc = _subInterfaces[i]->getConcVector();
// 			std::vector<double> subConc(_subDOFs[i], 0.0);
// 			for (auto pair : sparseConc[0][0][0]) {
// 				if (pair.first < _subDOFs[i])
// 					subConc[pair.first] = pair.second;
// 			}
// 			conc.push_back(subConc);
// 		}
// 		fullConc.push_back(conc);
// 	}
// 	// 1D
// 	else {
// 		for (auto j = 0; j < temperatures.size() - 2; j++) {
// 			// Loop on the grid points
// 			for (auto j = 0; j < temperatures.size() - 2; j++) {
// 				std::vector<std::vector<double>> conc;
// 				// Loop on the sub interfaces to get all the concentrations
// 				for (auto i = 0; i < _subInterfaces.size(); i++) {
// 					auto sparseConc = _subInterfaces[i]->getConcVector();
// 					std::vector<double> subConc(_subDOFs[i], 0.0);
// 					for (auto pair : sparseConc[0][0][j]) {
// 						if (pair.first < _subDOFs[i]) {
// 							subConc[pair.first] = pair.second;
// 						}
// 					}
// 					conc.push_back(subConc);
// 				}

// 				fullConc.push_back(conc);
// 			}
// 		}
// 	}

// 	// Print the result
// 	_networkInterface->outputData(
// 		_current_time, fullConc, std::max((int)temperatures.size() - 2, 1));
// }

// void
// XolotlNetworkProblem::externalSolve()
// {
// 	// Check that the next time is larger than the current one
// 	if (time() > _current_time) {
// 		double finalTime = 0.0, deltaTime = 0.0;
// 		if (dt() > _max_dt) {
// 			deltaTime = _max_dt;
// 			finalTime = _current_time + deltaTime;
// 		}
// 		else {
// 			deltaTime = dt();
// 			finalTime = time();
// 		}

// 		while (_current_time < time()) {
// 			// Transfer the temperature to the full network
// 			std::vector<double> temperatures;
// 			std::vector<double> depths;
// 			_subInterfaces[0]->getNetworkTemperature(temperatures, depths);
// 			std::vector<std::vector<std::vector<double>>> fullConc;

// 			// 0D
// 			if (temperatures.size() < 2) {
// 				std::vector<std::vector<double>> conc;
// 				// Loop on the sub interfaces to get all the concentrations
// 				for (auto i = 0; i < _subInterfaces.size(); i++) {
// 					auto sparseConc = _subInterfaces[i]->getConcVector();
// 					std::vector<double> subConc(_subDOFs[i], 0.0);
// 					for (auto pair : sparseConc[0][0][0]) {
// 						if (pair.first < _subDOFs[i]) {
// 							subConc[pair.first] = pair.second;
// 						}
// 					}
// 					conc.push_back(subConc);
// 				}

// 				fullConc.push_back(conc);

// 				// Compute the new rates
// 				std::vector<double> temperature = {temperatures[0]};
// 				std::vector<double> depth = {depths[0]};
// 				_networkInterface->setNetworkTemperature(temperature, depth);
// 				auto constantRates =
// 					_networkInterface->computeConstantRates(conc, 0);

// 				// Pass them
// 				for (auto i = 0; i < _subInterfaces.size(); i++) {
// 					_subInterfaces[i]->setConstantRates(constantRates[i], 0);
// 				}
// 			}
// 			// 1D
// 			else {
// 				// Loop on the grid points
// 				for (auto j = 0; j < temperatures.size() - 2; j++) {
// 					std::vector<std::vector<double>> conc;
// 					// Loop on the sub interfaces to get all the concentrations
// 					for (auto i = 0; i < _subInterfaces.size(); i++) {
// 						auto sparseConc = _subInterfaces[i]->getConcVector();
// 						std::vector<double> subConc(_subDOFs[i], 0.0);
// 						for (auto pair : sparseConc[0][0][j]) {
// 							if (pair.first < _subDOFs[i]) {
// 								subConc[pair.first] = pair.second;
// 							}
// 						}
// 						conc.push_back(subConc);
// 					}

// 					fullConc.push_back(conc);

// 					// Compute the new rates
// 					std::vector<double> temperature = {temperatures[j + 1]};
// 					std::vector<double> depth = {depths[j + 1]};
// 					_networkInterface->setNetworkTemperature(
// 						temperature, depth);
// 					auto constantRates =
// 						_networkInterface->computeConstantRates(conc, 0);

// 					// Pass them
// 					for (auto i = 0; i < _subInterfaces.size(); i++) {
// 						_subInterfaces[i]->setConstantRates(
// 							constantRates[i], j + 1);
// 					}
// 				}
// 			}

// 			// Print the result
// 			_networkInterface->outputData(_current_time, fullConc,
// 				std::max((int)temperatures.size() - 2, 1));

// 			// Solve
// 			for (auto i = 0; i < _subInterfaces.size(); i++) {
// 				// Set the time we want to reach
// 				_subInterfaces[i]->setTimes(finalTime, deltaTime);

// 				// Run the solver
// 				_subInterfaces[i]->solveXolotl();
// 			}
// 			// Save the current time
// 			_current_time += deltaTime;
// 			finalTime += deltaTime;
// 		}
// 	}
// }
} // namespace interface
} // namespace xolotl
