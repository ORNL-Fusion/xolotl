#include <cassert>
#include <fstream>
#include <limits>

using namespace std::string_literals;

#include <xolotl/options/ConfOptions.h>
#include <xolotl/options/InvalidOptionValue.h>
#include <xolotl/options/Options.h>
#ifdef Xolotl_ENABLE_JSON
#include <xolotl/options/JSONOptions.h>
#endif
#include <xolotl/util/Filesystem.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/Tokenizer.h>

namespace xolotl
{
namespace options
{
Options::Options() :
	petscArg(""),
	networkFilename(""),
	tempHandlerName(""),
	tempParam{},
	tempProfileFilename(""),
	tempGridPower(2.5),
	fluxFlag(false),
	fluxAmplitude(0.0),
	fluxTimeProfileFlag(false),
	perfHandlerName(""),
	perfOutputYAMLFlag(false),
	vizHandlerName(""),
	materialName(""),
	initialConcentration(""),
	interfaceLocation(-1000.0),
	dimensionNumber(1),
	gridTypeName(""),
	gridParam{},
	gridFilename(""),
	subnetworksFlag(false),
	initialTimeStep(0.0),
	maxTimeStep(0.0),
	timeStepGrowthFactor(0.0),
	startTime(0.0),
	endTime(0.0),
	numTimeSteps(0),
	gbList(""),
	groupingMin(std::numeric_limits<int>::max()),
	groupingWidthA(1),
	groupingWidthB(0),
	sputteringYield(0.0),
	useHDF5Flag(true),
	maxImpurity(8),
	maxD(0),
	maxT(0),
	maxV(20),
	maxPureV(20),
	maxI(6),
	leftBoundary(1),
	rightBoundary(1),
	bottomBoundary(1),
	topBoundary(1),
	frontBoundary(1),
	backBoundary(1),
	xBC("mirror"),
	heatLossPortion(-1.0),
	burstingDepth(10.0),
	burstingFactor(0.1),
	rngUseSeed(false),
	rngSeed(0),
	rngPrintSeed(false),
	zeta(0.73),
	density(10.162795276841),
	pulseTime(0.0),
	pulseProportion(0.0),
	latticeParameter(-1.0),
	impurityRadius(-1.0),
	biasFactor(1.15),
	hydrogenFactor(0.25),
	xenonDiffusivity(-1.0),
	fissionYield(0.25),
	heVRatio(4.0),
	migrationThreshold(std::numeric_limits<double>::infinity()),
	basalPortion(0.1),
	transitionSize(325),
	cascadeDose(-1.0),
	cascadeEfficiency(0.0)
{
}

Options::~Options()
{
}

void
Options::printAll(std::ostream& os) const
{
	os << "networkFilename: \"" << networkFilename << "\"\n";
	os << "tempHandlerName: \"" << tempHandlerName << "\"\n";
	os << "tempParam: {" << tempParam[0] << ", " << tempParam[1] << "}\n";
	os << "tempGridPower: " << tempGridPower << '\n';
	os << "tempProfileFilename: \"" << tempProfileFilename << "\"\n";
	os << "fluxFlag: " << std::boolalpha << fluxFlag << '\n';
	os << "fluxAmplitude: " << fluxAmplitude << '\n';
	os << "fluxTimeProfileFlag: " << std::boolalpha << fluxTimeProfileFlag
	   << '\n';
	os << "fluxTimeProfileFilePath: \"" << fluxTimeProfileFilePath << "\"\n";
	os << "perfHandlerName: \"" << perfHandlerName << "\"\n";
	os << "perfOutputYAMLFlag: " << std::boolalpha << perfOutputYAMLFlag
	   << '\n';
	os << "subnetworksFlag: " << std::boolalpha << subnetworksFlag << '\n';
	os << "initialTimeStep: " << initialTimeStep << '\n';
	os << "maxTimeStep: " << maxTimeStep << '\n';
	os << "timeStepGrowthFactor: " << timeStepGrowthFactor << '\n';
	os << "startTime: " << startTime << '\n';
	os << "endTime: " << endTime << '\n';
	os << "numTimeSteps: " << numTimeSteps << '\n';
	os << "vizHandlerName: \"" << vizHandlerName << "\"\n";
	os << "materialName: \"" << materialName << "\"\n";
	os << "initialConcentration: \"" << initialConcentration << "\"\n";
	os << "zeta: " << zeta << '\n';
	os << "interfaceLocation: " << interfaceLocation << '\n';
	os << "dimensionNumber: " << dimensionNumber << '\n';
	os << "gridTypeName: \"" << gridTypeName << "\"\n";
	os << "gridParam: {" << gridParam[0] << ", " << gridParam[1] << ", "
	   << gridParam[2] << ", " << gridParam[3] << ", " << gridParam[4] << ", "
	   << gridParam[5] << "}\n";
	os << "gridFilename: \"" << gridFilename << "\"\n";
	os << "gbList: \"" << gbList << "\"\n";
	os << "groupingMin: " << groupingMin << '\n';
	os << "groupingWidthA: " << groupingWidthA << '\n';
	os << "groupingWidthB: " << groupingWidthB << '\n';
	os << "sputteringYield: " << sputteringYield << '\n';
	os << "useHDF5Flag: " << std::boolalpha << useHDF5Flag << '\n';
	os << "networkParams:";
	for (auto&& p : networkParams) {
		os << ' ' << p;
	}
	os << '\n';
	os << "maxImpurity: " << maxImpurity << '\n';
	os << "maxD: " << maxD << '\n';
	os << "maxT: " << maxT << '\n';
	os << "maxV: " << maxV << '\n';
	os << "maxI: " << maxI << '\n';
	os << "maxPureV: " << maxPureV << '\n';

	os << "leftBoundary: " << leftBoundary << '\n';
	os << "rightBoundary: " << rightBoundary << '\n';
	os << "bottomBoundary: " << bottomBoundary << '\n';
	os << "topBoundary: " << topBoundary << '\n';
	os << "frontBoundary: " << frontBoundary << '\n';
	os << "backBoundary: " << backBoundary << '\n';

	os << "xBC: \"" << xBC << "\"\n";

	os << "heatLossPortion: " << heatLossPortion << '\n';
	os << "burstingDepth: " << burstingDepth << '\n';
	os << "burstingFactor: " << burstingFactor << '\n';
	os << "rngSeed: " << rngSeed << '\n';
	os << "rngUseSeed: " << std::boolalpha << rngUseSeed << '\n';
	os << "rngPrintSeed: " << std::boolalpha << rngPrintSeed << '\n';

	os << "radiusMinSizes:";
	for (auto&& s : radiusMinSizes) {
		os << ' ' << s;
	}
	os << '\n';

	os << "density: " << density << '\n';
	os << "pulseTime: " << pulseTime << '\n';
	os << "pulseProportion: " << pulseProportion << '\n';
	os << "latticeParameter: " << latticeParameter << '\n';
	os << "impurityRadius: " << impurityRadius << '\n';
	os << "biasFactor: " << biasFactor << '\n';
	os << "hydrogenFactor: " << hydrogenFactor << '\n';
	os << "xenonDiffusivity: " << xenonDiffusivity << '\n';
	os << "fissionYield: " << fissionYield << '\n';
	os << "heVRatio: " << heVRatio << '\n';
	os << "migrationThreshold: " << migrationThreshold << '\n';

	os << "fluxDepthProfileFilePath: \"" << fluxDepthProfileFilePath.string()
	   << "\"\n";

	os << "basalPortion: " << basalPortion << '\n';
	os << "transitionSize: " << transitionSize << '\n';
	os << "cascadeDose: " << cascadeDose << '\n';
	os << "cascadeEfficiency: " << cascadeEfficiency << '\n';

	os << "processMap:";
	for (auto&& p : processMap) {
		if (p.second) {
			os << " \"" << p.first << "\"";
		}
	}
	os << '\n';

	os << "petscArg: \"" << petscArg << "\"\n";
}

void
Options::setNetworkParameters(const std::vector<IdType>& params)
{
	networkParams = params;

	// Reset values
	maxD = 0;
	maxT = 0;
	maxV = 20;
	maxI = 6;
	maxPureV = 20;

	maxImpurity = params[0];

	if (params.size() > 1) {
		// Set the deuterium size
		maxD = params[1];
		// Set the tritium size
		maxT = params[2];
		// Set the vacancy size
		maxV = params[3];
		// Set the interstitial size
		maxI = params[4];
		if (params.size() > 5) {
			// Set the pure V size
			maxPureV = params[5];
		}
		else {
			maxPureV = maxV;
		}
	}
}

void
Options::setNetworkParameters(const std::string& paramStr)
{
	setNetworkParameters(util::Tokenizer<IdType>{paramStr}());
}

void
Options::addProcess(const std::string& processKey)
{
	// Initialize map if necessary
	if (processMap.empty()) {
		setProcesses("");
	}

	// Look for the key
	if (auto it = processMap.find(processKey); it != processMap.end()) {
		// Switch the value to true in the map
		processMap[processKey] = true;
	}
	else {
		throw InvalidOptionValue(
			"Options: The process name is not known: " + processKey);
	}
}

void
Options::setTempParam(const std::string& paramStr)
{
	auto tokens = util::Tokenizer<double>{paramStr}();
	if (tokens.size() > 2) {
		throw InvalidOptionValue(
			"Options: too many temperature parameters (expect 2 or less)");
	}
	for (std::size_t i = 0; i < tokens.size(); ++i) {
		tempParam[i] = tokens[i];
	}
}

void
Options::checkTempProfileFilename() const
{
	// Check that the profile file exists
	std::ifstream inFile(tempProfileFilename);
	if (!inFile) {
		throw std::runtime_error("Options: could not open file containing "
								 "temperature profile data. Aborting!");
	}
}

void
Options::checkFluxTimeProfileFilename() const
{
	// Check that the profile file exists
	std::ifstream inFile(fluxTimeProfileFilePath.c_str());
	if (!inFile) {
		throw std::runtime_error("Options: could not open file containing flux "
								 "profile data. Aborting!");
	}
}

void
Options::checkPerfHandlerName() const
{
	static const std::string perfHandlers[] = {"dummy", "os", "papi"};
	if (std::find(begin(perfHandlers), end(perfHandlers), perfHandlerName) ==
		end(perfHandlers)) {
		throw InvalidOptionValue(
			"Options: could not understand the performance handler type: " +
			perfHandlerName + ". Aborting!");
	}
}

void
Options::checkVizHandlerName() const
{
	// Determine the type of handlers we are being asked to use
	if (!(vizHandlerName == "std" || vizHandlerName == "dummy")) {
		throw InvalidOptionValue(
			"Options: unrecognized argument in the visualization option "
			"handler: " +
			vizHandlerName + ". Aborting!");
	}
}

void
Options::setGridParam(const std::string& paramStr)
{
	// Break the argument into tokens.
	auto tokens = util::Tokenizer<double>{paramStr}();
	if (tokens.size() > 6) {
		throw InvalidOptionValue(
			"Options: too many grid parameters (expect 6 or less)");
	}
	for (std::size_t i = 0; i < tokens.size(); ++i) {
		gridParam[i] = tokens[i];
	}
}

void
Options::checkGridFilename() const
{
	// Check that the file exists
	std::ifstream inFile(gridFilename);
	if (!inFile) {
		throw std::runtime_error(
			"Options: could not open file containing grid data. Aborting!");
	}
}

void
Options::setRadiusMinSizes(const std::string& paramStr)
{
	// Break the argument into tokens.
	auto tokens = util::Tokenizer<int>{paramStr}();

	// Set the values
	for (int i = 0; i < tokens.size(); i++) {
		radiusMinSizes.push_back(tokens[i]);
	}
}

void
Options::setBoundaries(const std::string& paramStr)
{
	// Break the argument into tokens.
	auto tokens = util::Tokenizer<int>{paramStr}();

	// Set the left boundary
	leftBoundary = tokens[0];
	// Set the right boundary
	rightBoundary = tokens[1];
	if (tokens.size() > 2)
		// Set the bottom boundary
		bottomBoundary = tokens[2];
	if (tokens.size() > 3)
		// Set the top boundary
		topBoundary = tokens[3];
	if (tokens.size() > 4)
		// Set the front boundary
		frontBoundary = tokens[4];
	if (tokens.size() > 5)
		// Set the back boundary
		backBoundary = tokens[5];
}

void
Options::processRNGParam(const std::string& paramStr)
{
	// Break the argument into tokens.
	auto tokens = util::Tokenizer<>{paramStr}();
	size_t currIdx = 0;

	// Determine whether we should print the seed value.
	bool shouldPrintSeed = false;
	if (tokens[currIdx] == "print") {
		shouldPrintSeed = true;
		++currIdx;
	}
	rngPrintSeed = shouldPrintSeed;

	if (currIdx < tokens.size()) {
		// Convert arg to an integer.
		char* ep = NULL;
		auto useed = strtoul(tokens[currIdx].c_str(), &ep, 10);
		if (ep != (tokens[currIdx].c_str() + tokens[currIdx].length())) {
			throw InvalidOptionValue(
				"Options: Invalid random number generator seed, must "
				"be a non-negative integer. Aborting!");
		}
		setRNGSeed(useed);
	}
}

void
Options::setProcesses(const std::string& processList)
{
	// Initialize the map of processes
	processMap["reaction"] = false;
	processMap["diff"] = false;
	processMap["advec"] = false;
	processMap["modifiedTM"] = false;
	processMap["movingSurface"] = false;
	processMap["bursting"] = false;
	processMap["attenuation"] = false;
	processMap["resolution"] = false;
	processMap["heterogeneous"] = false;
	processMap["sink"] = false;
	processMap["soret"] = false;
	processMap["constant"] = false;
	processMap["noSolve"] = false;

	// Break the argument into tokens.
	auto tokens = util::Tokenizer<>{processList}();

	// Loop on the tokens
	for (auto&& token : tokens) {
		addProcess(token);
	}
}

void
Options::setCouplingTimeStepParams(const std::string& paramString)
{
	// Set parameters from tokenized list
	auto params = util::Tokenizer<double>{paramString}();

	if (params.size() != 6) {
		throw InvalidOptionValue(
			"Options: Must provide six (6) values for time step "
			"parameters. Aborting!");
	}

	initialTimeStep = params[0];
	maxTimeStep = params[1];
	timeStepGrowthFactor = params[2];
	startTime = params[3];
	endTime = params[4];
	if (params[5] <= 0) {
		throw InvalidOptionValue(
			"Options: maxSteps must be a positive value. Aborting!");
	}
	numTimeSteps = static_cast<IdType>(params[5]);
}

void
Options::setPulseParams(const std::string& paramStr)
{
	// Break the argument into tokens.
	auto tokens = util::Tokenizer<double>{paramStr}();

	pulseTime = tokens[0];
	pulseProportion = tokens[1];
}

void
Options::setGroupingParams(const std::string& paramString)
{
	// Break the argument into tokens.
	auto tokens = util::Tokenizer<int>{paramString}();

	// Set grouping minimum size
	groupingMin = tokens[0];
	// Set the grouping width in the first direction
	groupingWidthA = tokens[1];
	// Set the grouping width in the second direction
	if (tokens.size() > 2) {
		groupingWidthB = tokens[2];
	}
}

void
Options::appendPetscArg(const std::string& arg)
{
	if (!petscArg.empty()) {
		petscArg += " ";
	}
	petscArg += arg;
}

std::shared_ptr<IOptions>
createOptions(int argc, const char* argv[])
{
	if (argc < 2) {
		throw std::runtime_error(
			"Options: parameter file name must not be empty");
	}

	auto filePath = fs::path(argv[1]);
	auto ext = filePath.extension();

#ifdef Xolotl_ENABLE_JSON
	if (ext == ".json") {
		return std::make_shared<JSONOptions>();
	}
#endif

	return std::make_shared<ConfOptions>();
}
} // end namespace options
} // end namespace xolotl
