#include <cassert>
#include <fstream>
#include <limits>

using namespace std::string_literals;

#include <boost/program_options.hpp>

#include <xolotl/options/ConfOptions.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/Tokenizer.h>

namespace bpo = boost::program_options;

namespace xolotl
{
namespace options
{
ConfOptions::ConfOptions() = default;

ConfOptions::ConfOptions(const ConfOptions& other) = default;

ConfOptions::~ConfOptions()
{
}

std::shared_ptr<IOptions>
ConfOptions::makeCopy() const
{
	return std::make_shared<ConfOptions>(*this);
}

void
ConfOptions::readParams(int argc, const char* argv[])
{
	// Check that a file name is given
	if (argc < 2) {
		throw bpo::error("Options: parameter file name must not be empty");
	}

	XOLOTL_LOG_WARN
		<< "WARNING: (ConfOptions) Parsing parameter files in config format is "
		   "deprecated and may not be up-to-date. "
		   "Consider converting your parameter file to JSON format.";

	// The name of the parameter file
	std::string param_file;

	// Parse the command line options.
	bpo::options_description desc("Command line options");
	desc.add_options()("help", "show this help message")("parameterFile",
		bpo::value<std::string>(&param_file),
		"When running Xolotl the name of the parameter file should immediately "
		"follow the executable: xolotl/xolotl param.txt")("restartFile",
		bpo::value<std::string>(&restartFile),
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
		"debug,extra,info,warning,error).")("restartFile",
		bpo::value<std::string>(&restartFile),
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
		"by the distance in nm, for instance: X 3.0 Z 2.5 Z 10.0 .")(
		"useSubnetworks", bpo::value<bool>(&subnetworksFlag),
		"Should we distribute network across subnetworks?")(
		"couplingTimeStepParams", bpo::value<std::string>(),
		"This option allows the user to define the parameters that control the "
		"multi-instance time-stepping. "
		"To do so, simply write the values in order "
		"initialDt maxDt growthFactor startTime endTime maxSteps.")("grouping",
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
		"migrationThreshold", bpo::value<double>(&migrationThreshold),
		"Set a limit on the migration energy above which the diffusion will be "
		"ignored.")("fluxDepthProfileFilePath",
		bpo::value<fs::path>(&fluxDepthProfileFilePath),
		"The path to the custom flux profile file; the default is an empty "
		"string that will use the default material associated flux handler.")(
		"reactionFilePath", bpo::value<fs::path>(&reactionFilePath),
		"The path to the reaction rates file; the default is an empty "
		"string.")("basalPortion",
		bpo::value<double>(&basalPortion)->default_value(0.1),
		"The value of the basal portion generated for each V (0.1 by "
		"default).")("transitionSize",
		bpo::value<int>(&transitionSize)->default_value(325),
		"The value for the transition within a type of cluster, for instance "
		"basal (325 by "
		"default).")("cascadeDose",
		bpo::value<double>(&cascadeDose)->default_value(-1.0),
		"The value of the dose at which the cascade overlap effect takes "
		"effect, if negative there won't be an effect (-1.0 by "
		"default).")("cascadeEfficiency",
		bpo::value<double>(&cascadeEfficiency)->default_value(0.0),
		"The value of the remaining efficiency once the overlap effect started "
		"(0.0 by "
		"default).");

	bpo::options_description visible("Allowed options");
	visible.add(desc).add(config);

	if (opts.count("help")) {
		std::cout << visible << std::endl;
		return;
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
		setTempParam(opts["tempParam"].as<std::string>());
	}

	if (opts.count("tempFile")) {
		checkTempProfileFilename();
	}

	// Take care of the flux
	if (opts.count("flux")) {
		fluxFlag = true;
	}
	if (opts.count("fluxFile")) {
		checkFluxTimeProfileFilename();
		// Set the flag to use a flux profile to true
		fluxTimeProfileFlag = true;
	}

	// Take care of the performance handler
	if (opts.count("perfHandler")) {
		checkPerfHandlerName();
	}

	// Take care of the visualization handler
	if (opts.count("vizHandler")) {
		checkVizHandlerName();
	}

	// Take care of the grid
	if (opts.count("gridParam")) {
		setGridParam(opts["gridParam"].as<std::string>());
	}

	if (opts.count("gridFile")) {
		checkGridFilename();
	}

	// Take care of the processes
	if (opts.count("process")) {
		setProcesses(opts["process"].as<std::string>());
	}

	// Take care of multi-instance params
	if (opts.count("couplingTimeStepParams")) {
		setCouplingTimeStepParams(
			opts["couplingTimeStepParams"].as<std::string>());
	}

	// Take care of the gouping
	if (opts.count("grouping")) {
		setGroupingParams(opts["grouping"].as<std::string>());
	}

	// Take care of the network parameters
	if (opts.count("netParam")) {
		// Set the flag to not use the HDF5 file
		useHDF5Flag = false;
		// Set parameters from tokenized list
		setNetworkParameters(opts["netParam"].as<std::string>());
	}

	// Take care of the radius minimum size
	if (opts.count("radiusSize")) {
		setRadiusMinSizes(opts["radiusSize"].as<std::string>());
	}

	// Take care of the boundary conditions
	if (opts.count("boundary")) {
		setBoundaries(opts["boundary"].as<std::string>());
	}

	// Take care of the rng
	if (opts.count("rng")) {
		processRNGParam(opts["rng"].as<std::string>());
	}

	// Take care of the flux pulse
	if (opts.count("pulse")) {
		setPulseParams(opts["pulse"].as<std::string>());
	}
}
} // namespace options
} // namespace xolotl
