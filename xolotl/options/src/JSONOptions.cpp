#include <fstream>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/algorithm/string/predicate.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <xolotl/options/JSONOptions.h>
#include <xolotl/options/detail/JSONElem.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>

using namespace std::string_literals;

namespace xolotl
{
namespace options
{
std::stringstream
stripComments(std::ifstream& ifs)
{
	std::stringstream ss;
	std::string line;
	bool cCommenting = false;
	while (std::getline(ifs, line)) {
		for (auto it = begin(line); it != end(line); ++it) {
			if (cCommenting) {
				if (*it == '*' && next(it) != end(line) && *next(it) == '/') {
					cCommenting = false;
					++it;
				}
				continue;
			}
			if (*it == '#') {
				break;
			}
			if (*it == '/') {
				if (next(it) != end(line) && *next(it) == '/') {
					break;
				}
				if (next(it) != end(line) && *next(it) == '*') {
					cCommenting = true;
					++it;
					continue;
				}
			}
			ss << *it;
		}
	}
	return ss;
}

template <typename TParam>
inline bool
checkSetParam(const boost::property_tree::iptree& tree,
	const std::string& label, TParam& param)
{
	auto wrap = tree.get_optional<TParam>(label);
	if (wrap) {
		param = wrap.get();
		return true;
	}
	return false;
}

template <typename T>
inline std::vector<T>
asVector(const boost::property_tree::iptree& node)
{
	std::vector<T> ret;
	for (auto&& item : node) {
		ret.push_back(item.second.get_value<T>());
	}
	return ret;
}

JSONOptions::JSONOptions() :
	_map(std::make_unique<boost::property_tree::iptree>())
{
}

JSONOptions::JSONOptions(const JSONOptions& other) :
	Options(other),
	_map(std::make_unique<boost::property_tree::iptree>(*other._map))
{
}

JSONOptions::~JSONOptions()
{
}

std::shared_ptr<IOptions>
JSONOptions::makeCopy() const
{
	return std::make_shared<JSONOptions>(*this);
}

detail::JSONElemVector
JSONOptions::defineHandlers()
{
#define JSON_ELEM_HANDLER [this, &tree](const std::string& name)

	using ElemType = detail::JSONElem::Type;

	auto& tree = *_map;

	detail::JSONElemVector handlers;

	handlers
		.add(
			"logLevel", ElemType::string,
			"Logging output threshold.\n"
			"{debug,extra,info,warning,error}\n"
			"(default = info)",
			JSON_ELEM_HANDLER {
				util::Log::setLevelThreshold(tree.get(name, "info"));
			})
		.add(
			"restartFile", ElemType::string,
			"The HDF5 file to use for restart.",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, restartFile); })
		.add(
			"tempHandler", ElemType::string,
			"Temperature handler to use.\n"
			"{constant,gradient,heat,ELM,profile}\n"
			"(default = constant)",
			JSON_ELEM_HANDLER { tempHandlerName = tree.get(name, "constant"); })
		.add(
			"tempParam", ElemType::list_real,
			"At most two parameters for temperature handler.\n"
			"Alternatives:\n"
			"constant -> temp;\n"
			"gradient -> surfaceTemp bulkTemp;\n"
			"heat -> heatFlux bulkTemp;\n"
			"ELM -> bulkTemp",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					auto node = tree.get_child(name);
					if (node.empty()) {
						setTempParam(node.get_value<std::string>());
					}
					else {
						setTempParam(asVector<double>(node));
					}
				}
			})
		.add(
			"tempFile", ElemType::string,
			"A temperature profile is given by the specified file, "
			"then linear interpolation is used to fit the data.\n"
			"NOTE: no need for tempParam here.",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					tempProfileFilename = tree.get<std::string>(name);
					checkTempProfileFilename();
				}
			})
		.add(
			"tempGridPower", ElemType::real,
			"The value of the power to use to create the temperature grid "
			"spacing, only used if heat temperature handler is used.\n"
			"(default = 2.5)",
			JSON_ELEM_HANDLER {
				checkSetParam(tree, "tempGridPower", tempGridPower);
			})
		.add(
			"flux", ElemType::real,
			"The value of the incoming flux in #/nm2/s.\n"
			"If the Fuel case is used it actually corresponds to the fission "
			"rate in #/nm3/s.",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					fluxFlag = true;
					fluxAmplitude = tree.get<double>(name);
				}
			})
		.add(
			"fluxFile", ElemType::string,
			"A time profile for the flux is given by the specified file, "
			"then linear interpolation is used to fit the data.\n"
			"NOTE: If a flux profile file is given, "
			"a constant flux should NOT be given",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					fluxTimeProfileFilePath = tree.get<std::string>(name);
					checkFluxTimeProfileFilename();
					fluxTimeProfileFlag = true;
				}
			})
		.add(
			"perfHandler", ElemType::string,
			"Which set of performance handlers to use.\n"
			"{dummy,os,papi}\n"
			"(default = os)",
			JSON_ELEM_HANDLER {
				perfHandlerName = tree.get(name, "os");
				checkPerfHandlerName();
			})
		.add(
			"perfOutputYAML", ElemType::boolean,
			"Should we write the performance report to a YAML file?",
			JSON_ELEM_HANDLER {
				checkSetParam(tree, name, perfOutputYAMLFlag);
			})
		.add(
			"vizHandler", ElemType::string,
			"Which set of handlers to use for the visualization.\n"
			"{std,dummy}\n"
			"(default = dummy)",
			JSON_ELEM_HANDLER {
				vizHandlerName = tree.get(name, "dummy");
				checkVizHandlerName();
			})
		.add(
			"dimensions", ElemType::integer,
			"Number of dimensions for the simulation.",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, dimensionNumber); })
		.add(
			"material", ElemType::string,
			"Label for material being used.\n"
			"{W100, W110, W111, W211, Pulsed, Fuel, Fe, 800H, AlphaZr}",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, materialName); })
		.add(
			"initialConc", ElemType::string,
			"The name, size, and value of the initial concentration in the "
			"material.",
			JSON_ELEM_HANDLER {
				checkSetParam(tree, name, initialConcentration);
			})
		.add(
			"zeta", ElemType::real,
			"The value of the electronic stopping power in the material.\n"
			"(default = 0.73)",
			JSON_ELEM_HANDLER { checkSetParam(tree, "zeta", zeta); })
		.add(
			"interfaceLoc", ElemType::real,
			"The value (in nm) of the interface location between two "
			"materials.\n"
			"(default = -1000.0)",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, interfaceLocation); })
		.add(
			"gridType", ElemType::string,
			"Grid type to use along X.\n"
			"{uniform,nonuniform,geometric,cheby,read}\n"
			"(default = uniform)",
			JSON_ELEM_HANDLER { gridTypeName = tree.get(name, "uniform"); })
		.add(
			"gridParam", ElemType::list_real,
			"At most six parameters for the grid.\n"
			"Alternatives:\n"
			"uniform -> nX hX;\n"
			"nonuniform -> nX;\n"
			"geometric -> nX ratio\n"
			"cheby -> nX width\n"
			"The four additional parameters are for a uniform grid in Y and Z\n"
			"-> nY hY nZ hZ.",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					auto node = tree.get_child(name);
					if (node.empty()) {
						setGridParam(node.get_value<std::string>());
					}
					else {
						setGridParam(asVector<double>(node));
					}
				}
			})
		.add(
			"gridFile", ElemType::string,
			"A grid spacing is given by the specified file.\n"
			"NOTE: you need gridParam here only if in 2D/3D.",
			JSON_ELEM_HANDLER {
				if (checkSetParam(tree, name, gridFilename)) {
					checkGridFilename();
				}
			})
		.add(
			"process", ElemType::list_string,
			"List of all the processes to use in the simulation.\n"
			"{reaction, diff, advec, modifiedTM, movingSurface, "
			"bursting, attenuation, resolution, heterogeneous, "
			"sink, soret, constant, noSolve}",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					auto node = tree.get_child(name);
					if (node.empty()) {
						setProcesses(node.get_value<std::string>());
					}
					else {
						setProcesses(asVector<std::string>(node));
					}
				}
			})
		.add(
			"grain", ElemType::string,
			"This option allows the user to add GB in the X, Y, or Z "
			"directions. To do so, simply write the direction followed "
			"by the distance in nm.\n"
			"For instance: X 3.0 Z 2.5 Z 10.0",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, gbList); })
		.add(
			"useSubnetworks", ElemType::boolean,
			"Should we distribute network across subnetworks?",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, subnetworksFlag); })
		.add(
			"couplingTimeStepParams", ElemType::list_real,
			"This option allows the user to define the parameters that control "
			"the multi-instance time-stepping. "
			"To do so, simply write the values in order:\n"
			"initialDt maxDt growthFactor startTime endTime maxSteps",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					auto node = tree.get_child(name);
					if (node.empty()) {
						setCouplingTimeStepParams(
							node.get_value<std::string>());
					}
					else {
						setCouplingTimeStepParams(asVector<double>(node));
					}
				}
			})
		.add(
			"grouping", ElemType::list_integer,
			"The grouping parameters: the first integer is the size at which "
			"the grouping starts (HeV clusters in the PSI case, Xe in the NE "
			"case), the second is the first width of the groups (He for PSI, "
			"Xe for NE), and the third one in the second width of the groups "
			"(V for PSI).",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					auto node = tree.get_child(name);
					if (node.empty()) {
						setGroupingParams(node.get_value<std::string>());
					}
					else {
						setGroupingParams(asVector<int>(node));
					}
				}
			})
		.add(
			"sputtering", ElemType::real,
			"The sputtering yield (in atoms/ion) that will be used.",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, sputteringYield); })
		.add(
			"netParam", ElemType::list_integer,
			"This option allows the user to define the boundaries of the "
			"network. To do so, simply write the values in order:\n"
			"maxHe/Xe/Basal maxD maxT maxV maxI maxPureV",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					// Set the flag to not use the HDF5 file
					useHDF5Flag = false;
					auto node = tree.get_child(name);
					if (node.empty()) {
						setNetworkParameters(node.get_value<std::string>());
					}
					else {
						setNetworkParameters(asVector<IdType>(node));
					}
				}
			})
		.add(
			"radiusSize", ElemType::list_integer,
			"This option allows the user to set a minimum size for the "
			"computation for the average radii, in the same order as the "
			"netParam option (default is 0).",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					auto node = tree.get_child(name);
					if (node.empty()) {
						setRadiusMinSizes(node.get_value<std::string>());
					}
					else {
						setRadiusMinSizes(asVector<int>(node));
					}
				}
			})
		.add(
			"boundary", ElemType::list_integer,
			"This option allows the user to choose the boundary conditions. "
			"The first one corresponds to the left side (surface) and the "
			"second one to the right (bulk), then two for Y and two for Z.\n"
			"0 means mirror or periodic, 1 means free surface.",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					auto node = tree.get_child(name);
					if (node.empty()) {
						setBoundaries(node.get_value<std::string>());
					}
					else {
						setBoundaries(asVector<int>(node));
					}
				}
			})
		.add(
			"xBCType", ElemType::string,
			"The boundary conditions to use in the X direction.\n"
			"{mirror,periodic,robin (for temperature)}\n"
			"(default = mirror)",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, xBC); })
		.add(
			"heatLossPortion", ElemType::real,
			"The portion of heat lost in the bulk.\n"
			"(default = -1.0)",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, heatLossPortion); })
		.add(
			"burstingDepth", ElemType::real,
			"The depth (in nm) after which there is an exponential decrease in "
			"the probability of bursting.\n"
			"(default = 10.0)",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, burstingDepth); })
		.add(
			"burstingFactor", ElemType::real,
			"This option allows the user to set the factor used in computing "
			"the likelihood of a bursting event.",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, burstingFactor); })
		.add(
			"rng", ElemType::string,
			"Allows user to specify seed used to initialize random number "
			"generator (default = determined from current time) and "
			"whether each process should print the seed value "
			"it uses (default = don't print).",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					processRNGParam(tree.get<std::string>(name));
				}
			})
		.add(
			"density", ElemType::real,
			"Sets a density in nm-3 for the number of xenon per volume in a "
			"bubble for the NE case.\n"
			"(default = 10.162795276841)",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, density); })
		.add(
			"pulse", ElemType::string,
			"The total length of the pulse (in s) if the Pulsed material is "
			"used, and the proportion of it that is ON.",
			JSON_ELEM_HANDLER {
				if (tree.count(name)) {
					setPulseParams(tree.get<std::string>(name));
				}
			})
		.add(
			"lattice", ElemType::real, "The length of the lattice side in nm.",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, latticeParameter); })
		.add(
			"impurityRadius", ElemType::real,
			"The radius of the main impurity (He or Xe) in nm.",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, impurityRadius); })
		.add(
			"biasFactor", ElemType::real,
			"This option allows the user to set the bias factor reflecting the "
			"fact that interstitial clusters have a larger surrounding strain "
			"field.",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, biasFactor); })
		.add(
			"hydrogenFactor", ElemType::real,
			"The factor between the size of He and H.",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, hydrogenFactor); })
		.add(
			"xenonDiffusivity", ElemType::real,
			"The diffusion coefficient for xenon in nm2 s-1.",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, xenonDiffusivity); })
		.add(
			"fissionYield", ElemType::real,
			"The number of xenon created for each fission.\n"
			"(default = 0.25)",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, fissionYield); })
		.add(
			"heVRatio", ElemType::real,
			"The number of He atoms allowed per V in a bubble.",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, heVRatio); })
		.add(
			"migrationThreshold", ElemType::real,
			"Set a limit on the migration energy above which the diffusion "
			"will be ignored.",
			JSON_ELEM_HANDLER {
				checkSetParam(tree, name, migrationThreshold);
			})
		.add(
			"fluxDepthProfileFilePath", ElemType::string,
			"The path to the custom flux profile file; the default is an empty "
			"string that will use the default material associated flux "
			"handler.",
			JSON_ELEM_HANDLER {
				checkSetParam(tree, name, fluxDepthProfileFilePath);
			})
		.add(
			"reactionFilePath", ElemType::string,
			"The path to the reaction rates file.\n"
			"(default = \"\")",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, reactionFilePath); })
		.add(
			"basalPortion", ElemType::real,
			"The value of the basal portion generated for each V.\n"
			"(default = 0.1)",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, basalPortion); })
		.add(
			"transitionSize", ElemType::integer,
			"The value for the transition within a type of cluster, for "
			"instance basal.\n"
			"(default = 325)",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, transitionSize); })
		.add(
			"cascadeDose", ElemType::real,
			"The value of the dose at which the cascade overlap effect takes "
			"effect, if negative there won't be an effect.\n"
			"(default = -1.0)",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, cascadeDose); })
		.add(
			"cascadeEfficiency", ElemType::real,
			"The value of the remaining efficiency once the overlap effect "
			"started.\n"
			"(default = 0.0)",
			JSON_ELEM_HANDLER { checkSetParam(tree, name, cascadeEfficiency); })
		.add(
			"petscArgs", ElemType::list_string,
			"List of strings to be joined verbatim and passed to PETSc.",
			JSON_ELEM_HANDLER {
				if (tree.count("petscArgs")) {
					for (auto&& elem : tree.get_child("petscArgs")) {
						appendPetscArg(elem.second.data());
					}
				}
			});

	return handlers;

#undef JSON_ELEM_HANDLER
}

void
JSONOptions::printHelp(std::ostream& os)
{
	defineHandlers().print(std::cout);
	std::cout
		<< "NOTE: a LIST parameter may be specified either as a JSON list\n"
		   "      or as a space separated string list.\n\n";
}

void
JSONOptions::readParams(int argc, const char* argv[])
{
	if (argc < 2) {
		throw std::runtime_error(
			"Options: parameter file name must not be empty");
	}

	std::ifstream ifs(argv[1]);
	if (!ifs) {
		throw std::runtime_error("Unable to open file: "s + argv[1]);
	}

	auto ss = stripComments(ifs);
	boost::property_tree::read_json(ss, *_map);

	defineHandlers().processParams();
}
} // namespace options
} // namespace xolotl
