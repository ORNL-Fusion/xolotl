#include <cassert>
#include <fstream>
#include <limits>

using namespace std::string_literals;

#include <boost/program_options.hpp>

#include <xolotl/options/Options.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/Tokenizer.h>

namespace bpo = boost::program_options;

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
	cascadeDose(-1.0),
	cascadeEfficiency(0.0),
	barrierEnergy(1.0)
{
	return;
}

Options::~Options(void)
{
}

void
Options::readParams(int argc, const char* argv[])
{
	// Check that a file name is given
	if (argc < 2) {
		throw bpo::error("Options: parameter file name must not be empty");
	}

	// The name of the parameter file
	std::string param_file;

	// Parse the command line options.
	bpo::options_description desc("Command line options");
	desc.add_options()("help", "show this help message")("parameterFile",
		bpo::value<std::string>(&param_file),
		"When running Xolotl the name of the parameter file should immediately "
		"follow the executable: xolotl/xolotl param.txt")("networkFile",
		bpo::value<std::string>(&networkFilename),
		"The HDF5 file to use for restart.");

	bpo::positional_options_description p;
	p.add("parameterFile", 1);

	bpo::variables_map opts;

	bpo::store(
		bpo::command_line_parser(argc, argv).options(desc).positional(p).run(),
		opts);
	bpo::notify(opts);

	// Declare a group of options that will be
	// allowed in config file
	bpo::options_description config("Parameters");
	config.add_options()("logLevel",
		bpo::value<std::string>()->default_value("info"),
		"Logging output threshold. (default = info; available "
		"debug,extra,info,warning,error).")("networkFile",
		bpo::value<std::string>(&networkFilename),
		"The HDF5 file to use for restart.")("tempHandler",
		bpo::value<std::string>(&tempHandlerName)->default_value("constant"),
		"Temperature handler to use. (default = constant; available "
		"constant,gradient,heat,ELM,profile")("tempParam",
		bpo::value<std::string>(),
		"At most two parameters for temperature handler. Alternatives:"
		"constant -> temp; "
		"gradient -> surfaceTemp bulkTemp; "
		"heat -> heatFlux bulkTemp; "
		"ELM -> bulkTemp")("tempFile",
		bpo::value<std::string>(&tempProfileFilename),
		"A temperature profile is given by the specified file, "
		"then linear interpolation is used to fit the data."
		" NOTE: no need for tempParam here.")("tempGridPower",
		bpo::value<double>(&tempGridPower),
		"The value of the power to use to create the temperature grid spacing, "
		"only used if heat temperature handler is used. (default = 2.5).")(
		"flux", bpo::value<double>(&fluxAmplitude),
		"The value of the incoming flux in #/nm2/s. If the Fuel case is used "
		"it actually corresponds to the fission rate in #/nm3/s.")("fluxFile",
		bpo::value<std::string>(&fluxTimeProfileFilePath),
		"A time profile for the flux is given by the specified file, "
		"then linear interpolation is used to fit the data."
		"(NOTE: If a flux profile file is given, "
		"a constant flux should NOT be given)")("perfHandler",
		bpo::value<std::string>(&perfHandlerName)->default_value("os"),
		"Which set of performance handlers to use. (default = os, available "
		"dummy,os,papi).")("perfOutputYAML",
		bpo::value<bool>(&perfOutputYAMLFlag),
		"Should we write the performance report to a YAML file?")("vizHandler",
		bpo::value<std::string>(&vizHandlerName)->default_value("dummy"),
		"Which set of handlers to use for the visualization. (default = dummy, "
		"available std,dummy).")("dimensions",
		bpo::value<int>(&dimensionNumber),
		"Number of dimensions for the simulation.")("material",
		bpo::value<std::string>(&materialName),
		"The material options are as follows: {W100, W110, W111, "
		"W211, Pulsed, Fuel, Fe, 800H, AlphaZr}.")("initialConc",
		bpo::value<std::string>(&initialConcentration),
		"The name, size, and value of the initial concentration in the "
		"material.")("zeta", bpo::value<double>(&zeta)->default_value(0.73),
		"The value of the electronic stopping power in the material (0.73 by "
		"default).")("interfaceLoc", bpo::value<double>(&interfaceLocation),
		"The value (in nm) of the interface location between two materials "
		"(-1000.0 nm by default).")("gridType",
		bpo::value<std::string>(&gridTypeName)->default_value("uniform"),
		"Grid type to use along X. (default = uniform; available "
		"uniform,nonuniform,geometric,cheby,read")("gridParam",
		bpo::value<std::string>(),
		"At most six parameters for the grid. Alternatives:"
		"uniform -> nX hX; "
		"nonuniform -> nX; "
		"geometric -> nX ratio"
		"cheby -> nX width"
		"The four additional parameters are for a uniform grid in Y and Z -> "
		"nY hY nZ hZ.")("gridFile", bpo::value<std::string>(&gridFilename),
		"A grid spacing is given by the specified file."
		" NOTE: you need gridParam here only if in 2D/3D.")("petscArgs",
		bpo::value<std::string>(&petscArg),
		"All the arguments that will be given to PETSc.")("process",
		bpo::value<std::string>(),
		"List of all the processes to use in the simulation (reaction, diff, "
		"advec, modifiedTM, movingSurface, bursting, attenuation, resolution, "
		"heterogeneous, sink, soret, constant, noSolve).")("grain",
		bpo::value<std::string>(&gbList),
		"This option allows the user to add GB in the X, Y, or Z directions. "
		"To do so, simply write the direction followed "
		"by the distance in nm, for instance: X 3.0 Z 2.5 Z 10.0 .")("grouping",
		bpo::value<std::string>(),
		"The grouping parameters: the first integer is the size at which the "
		"grouping starts (HeV clusters in the PSI case, Xe in the NE case), "
		"the second is the first width of the groups (He for PSI, Xe for NE), "
		"and the third one in the second width of the groups (V for PSI).")(
		"sputtering", bpo::value<double>(&sputteringYield),
		"The sputtering yield (in atoms/ion) that will be used.")("netParam",
		bpo::value<std::string>(),
		"This option allows the user to define the boundaries of the network. "
		"To do so, simply write the values in order "
		"maxHe/Xe/Basal maxD maxT maxV maxI maxPureV.")("radiusSize",
		bpo::value<std::string>(),
		"This option allows the user to set a minimum size for the computation "
		"for the average radii, in the same order as the netParam option "
		"(default is 0).")("boundary", bpo::value<std::string>(),
		"This option allows the user to choose the boundary conditions. "
		"The first one correspond to the left side (surface) "
		"and second one to the right (bulk), "
		"then two for Y and two for Z. "
		"0 means mirror or periodic, 1 means free surface.")("xBCType",
		bpo::value<std::string>(&xBC),
		"The boundary conditions to use in the X direction, mirror (default), "
		"periodic, or robin (for temperature).")("heatLossPortion",
		bpo::value<double>(&heatLossPortion),
		"The portion of heat lost in the bulk (-1.0 by default).")(
		"burstingDepth", bpo::value<double>(&burstingDepth),
		"The depth (in nm) after which there is an exponential decrease in the "
		"probability of bursting (10.0 nm if nothing is specified).")(
		"burstingFactor", bpo::value<double>(&burstingFactor),
		"This option allows the user to set the factor used in computing the "
		"likelihood of a bursting event.")("rng", bpo::value<std::string>(),
		"Allows user to specify seed used to initialize random number "
		"generator (default = determined from current time) and "
		"whether each process should print the seed value "
		"it uses (default = don't print).")("density",
		bpo::value<double>(&density),
		"Sets a density in nm-3 for the number of xenon per volume in a bubble "
		"for the NE case (default is 10.162795276841 nm-3 as before).")("pulse",
		bpo::value<std::string>(),
		"The total length of the pulse (in s) if the Pulsed material is used, "
		"and the proportion of it that is "
		"ON.")("lattice", bpo::value<double>(&latticeParameter),
		"The length of the lattice side in nm.")("impurityRadius",
		bpo::value<double>(&impurityRadius),
		"The radius of the main impurity (He or Xe) in nm.")("biasFactor",
		bpo::value<double>(&biasFactor),
		"This option allows the user to set the bias factor reflecting the "
		"fact that interstitial "
		"clusters have a larger surrounding strain field.")("hydrogenFactor",
		bpo::value<double>(&hydrogenFactor),
		"The factor between the size of He and H.")("xenonDiffusivity",
		bpo::value<double>(&xenonDiffusivity),
		"The diffusion coefficient for xenon in nm2 s-1.")("fissionYield",
		bpo::value<double>(&fissionYield),
		"The number of xenon created for each fission (default is 0.25).")(
		"heVRatio", bpo::value<double>(&heVRatio),
		"The number of He atoms allowed per V in a bubble.")(
		"migrationThreshold", bpo::value<double>(&migrationThreshold),
		"Set a limit on the migration energy above which the diffusion will be "
		"ignored.")("fluxDepthProfileFilePath",
		bpo::value<fs::path>(&fluxDepthProfileFilePath),
		"The path to the custom flux profile file; the default is an empty "
		"string that will use the default material associated flux handler.")(
		"basalPortion", bpo::value<double>(&basalPortion)->default_value(0.1),
		"The value of the basal portion generated for each V (0.1 by "
		"default).")("transitionSize", bpo::value<std::string>(),
		"The size at which there is a transition, for basal for instance or "
		"vacancy and interstitial transformation for alloy (it can take "
		"multiple values.")("cascadeDose",
		bpo::value<double>(&cascadeDose)->default_value(-1.0),
		"The value of the dose at which the cascade overlap effect takes "
		"effect, if negative there won't be an effect (-1.0 by "
		"default).")("cascadeEfficiency",
		bpo::value<double>(&cascadeEfficiency)->default_value(0.0),
		"The value of the remaining efficiency once the overlap effect started "
		"(0.0 by "
		"default).")("barrierEnergy",
		bpo::value<double>(&barrierEnergy)->default_value(1.0),
		"The value of the energy barrier to form 100 loops "
		"(1.0 eV by "
		"default).");

	bpo::options_description visible("Allowed options");
	visible.add(desc).add(config);

	if (opts.count("help")) {
		std::cout << visible << std::endl;
		std::exit(EXIT_SUCCESS);
	}

	// Check that the file exist
	if (!std::ifstream{argv[1]}) {
		std::cout << visible << std::endl;
		throw bpo::error("Options: unable to open parameter file: "s + argv[1]);
	}

	std::ifstream ifs(param_file);
	if (!ifs) {
		throw bpo::error(
			"Options: unable to open parameter file: " + param_file);
	}
	store(parse_config_file(ifs, config), opts);
	notify(opts);

	util::Log::setLevelThreshold(opts["logLevel"].as<std::string>());

	// Take care of the temperature
	if (opts.count("tempParam")) {
		// Break the argument into tokens.
		auto tokens =
			util::Tokenizer<double>{opts["tempParam"].as<std::string>()}();
		if (tokens.size() > 2) {
			throw bpo::error(
				"Options: too many temperature parameters (expect 2 or less)");
		}
		for (std::size_t i = 0; i < tokens.size(); ++i) {
			tempParam[i] = tokens[i];
		}
	}

	if (opts.count("tempFile")) {
		// Check that the profile file exists
		std::ifstream inFile(tempProfileFilename);
		if (!inFile) {
			throw bpo::error("Options: could not open file containing "
							 "temperature profile data. Aborting!");
		}
	}

	// Take care of the flux
	if (opts.count("flux")) {
		fluxFlag = true;
	}
	if (opts.count("fluxFile")) {
		// Check that the profile file exists
		std::ifstream inFile(fluxTimeProfileFilePath.c_str());
		if (!inFile) {
			throw bpo::error("Options: could not open file containing flux "
							 "profile data. Aborting!");
		}
		else {
			// Set the flag to use a flux profile to true
			fluxTimeProfileFlag = true;
		}
	}

	// Take care of the performance handler
	if (opts.count("perfHandler")) {
		std::string perfHandlers[] = {"dummy", "os", "papi"};
		if (std::find(begin(perfHandlers), end(perfHandlers),
				perfHandlerName) == end(perfHandlers)) {
			throw bpo::invalid_option_value(
				"Options: could not understand the performance handler type. "
				"Aborting!");
		}
	}

	// Take care of the visualization handler
	if (opts.count("vizHandler")) {
		// Determine the type of handlers we are being asked to use
		if (!(vizHandlerName == "std" || vizHandlerName == "dummy")) {
			throw bpo::invalid_option_value(
				"Options: unrecognized argument in the visualization option "
				"handler. Aborting!");
		}
	}

	// Take care of the grid
	if (opts.count("gridParam")) {
		// Break the argument into tokens.
		auto tokens =
			util::Tokenizer<double>{opts["gridParam"].as<std::string>()}();
		if (tokens.size() > 6) {
			throw bpo::invalid_option_value(
				"Options: too many grid parameters (expect 6 or less)");
		}
		for (std::size_t i = 0; i < tokens.size(); ++i) {
			gridParam[i] = tokens[i];
		}
	}

	if (opts.count("gridFile")) {
		// Check that the file exists
		std::ifstream inFile(gridFilename);
		if (!inFile) {
			throw bpo::invalid_option_value(
				"Options: could not open file containing grid data. Aborting!");
		}
	}

	// Take care of the radius minimum size
	if (opts.count("radiusSize")) {
		// Break the argument into tokens.
		auto tokens =
			util::Tokenizer<int>{opts["radiusSize"].as<std::string>()}();

		// Set the values
		for (int i = 0; i < tokens.size(); i++) {
			radiusMinSizes.push_back(tokens[i]);
		}
	}

	// Take care of the processes
	if (opts.count("process")) {
		// Break the argument into tokens.
		auto tokens = util::Tokenizer<>{opts["process"].as<std::string>()}();

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

		// Loop on the tokens
		for (int i = 0; i < tokens.size(); ++i) {
			// Look for the key
			if (processMap.find(tokens[i]) == processMap.end()) {
				throw bpo::invalid_option_value(
					"Options: The process name is not known: "s + tokens[i] +
					". Aborting!");
			}
			else {
				// Switch the value to true in the map
				processMap[tokens[i]] = true;
			}
		}
	}

	// Take care of the gouping
	if (opts.count("grouping")) {
		// Break the argument into tokens.
		auto tokens =
			util::Tokenizer<int>{opts["grouping"].as<std::string>()}();

		// Set grouping minimum size
		groupingMin = tokens[0];
		// Set the grouping width in the first direction
		groupingWidthA = tokens[1];
		// Set the grouping width in the second direction
		if (tokens.size() > 2)
			groupingWidthB = tokens[2];
	}

	// Take care of the network parameters
	if (opts.count("netParam")) {
		// Break the argument into tokens.
		auto tokens = util::Tokenizer<>{opts["netParam"].as<std::string>()}();

		// Set the flag to not use the HDF5 file
		useHDF5Flag = false;

		// Set the value for the impurities
		maxImpurity = strtol(tokens[0].c_str(), NULL, 10);

		// Check if we have other values
		if (tokens.size() > 1) {
			// Set the deuterium size
			maxD = strtol(tokens[1].c_str(), NULL, 10);
			// Set the tritium size
			maxT = strtol(tokens[2].c_str(), NULL, 10);
			// Set the vacancy size
			maxV = strtol(tokens[3].c_str(), NULL, 10);
			// Set the interstitial size
			maxI = strtol(tokens[4].c_str(), NULL, 10);
			if (tokens.size() > 5) {
				// Set the pure V size
				maxPureV = strtol(tokens[5].c_str(), NULL, 10);
			}
			else
				maxPureV = maxV;
		}
	}

	// Take care of the boundary conditions
	if (opts.count("boundary")) {
		// Break the argument into tokens.
		auto tokens =
			util::Tokenizer<int>{opts["boundary"].as<std::string>()}();

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

	// Take care of the rng
	if (opts.count("rng")) {
		// Break the argument into tokens.
		auto tokens = util::Tokenizer<>{opts["rng"].as<std::string>()}();
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
				throw bpo::invalid_option_value(
					"Options: Invalid random number generator seed, must "
					"be a non-negative integer. Aborting!");
			}
			setRNGSeed(useed);
		}
	}

	// Take care of the flux pulse
	if (opts.count("pulse")) {
		// Break the argument into tokens.
		auto tokens =
			util::Tokenizer<double>{opts["pulse"].as<std::string>()}();

		pulseTime = tokens[0];
		pulseProportion = tokens[1];
	}

	// Take care of the transition size
	if (opts.count("transitionSize")) {
		// Break the argument into tokens.
		auto tokens =
			util::Tokenizer<int>{opts["transitionSize"].as<std::string>()}();
		for (std::size_t i = 0; i < tokens.size(); ++i) {
			transitionSize.push_back(tokens[i]);
		}
	}
	else
		transitionSize.push_back(325);
}

} // end namespace options
} // end namespace xolotl
