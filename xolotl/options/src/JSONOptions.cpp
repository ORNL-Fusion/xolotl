#include <fstream>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/algorithm/string/predicate.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <xolotl/options/JSONOptions.h>
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

JSONOptions::JSONOptions() = default;

JSONOptions::JSONOptions(const JSONOptions& other) :
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

	_map = std::make_unique<boost::property_tree::iptree>();
	auto& tree = *_map;

	auto ss = stripComments(ifs);
	boost::property_tree::read_json(ss, tree);

	util::Log::setLevelThreshold(tree.get("logLevel", "info"));

	checkSetParam(tree, "networkFile", networkFilename);

	tempHandlerName = tree.get("tempHandler", "constant");

	if (tree.count("tempParam")) {
        auto node = tree.get_child("tempParam");
        if (node.empty()) {
            setTempParam(node.get_value<std::string>());
        }
        else {
            setTempParam(asVector<double>(node));
        }
	}

	if (tree.count("tempFile")) {
		tempProfileFilename = tree.get<std::string>("tempFile");
		checkTempProfileFilename();
	}

	checkSetParam(tree, "tempGridPower", tempGridPower);

	if (tree.count("flux")) {
		fluxFlag = true;
		fluxAmplitude = tree.get<double>("flux");
	}
	if (tree.count("fluxFile")) {
		fluxTimeProfileFilePath = tree.get<std::string>("fluxFile");
		checkFluxTimeProfileFilename();
		fluxTimeProfileFlag = true;
	}

	perfHandlerName = tree.get("perfHandler", "os");
	checkPerfHandlerName();

	checkSetParam(tree, "perfOutputYAML", perfOutputYAMLFlag);

	vizHandlerName = tree.get("vizHandler", "dummy");
	checkVizHandlerName();

	checkSetParam(tree, "dimensions", dimensionNumber);

	checkSetParam(tree, "material", materialName);

	checkSetParam(tree, "initialConc", initialConcentration);

	checkSetParam(tree, "zeta", zeta);

	checkSetParam(tree, "interfaceLoc", interfaceLocation);

	gridTypeName = tree.get("gridType", "uniform");

	if (tree.count("gridParam")) {
		auto node = tree.get_child("gridParam");
		if (node.empty()) {
			setGridParam(node.get_value<std::string>());
		}
		else {
			setGridParam(asVector<double>(node));
		}
	}

	if (checkSetParam(tree, "gridFile", gridFilename)) {
		checkGridFilename();
	}

	if (tree.count("process")) {
		auto node = tree.get_child("process");
		if (node.empty()) {
			setProcesses(node.get_value<std::string>());
		}
		else {
			setProcesses(asVector<std::string>(node));
		}
	}

	checkSetParam(tree, "grain", gbList);

	checkSetParam(tree, "useSubnetworks", subnetworksFlag);

	if (tree.count("couplingTimeStepParams")) {
		auto node = tree.get_child("couplingTimeStepParams");
		if (node.empty()) {
			setCouplingTimeStepParams(node.get_value<std::string>());
		}
		else {
			setCouplingTimeStepParams(asVector<double>(node));
		}
	}

	if (tree.count("grouping")) {
		auto node = tree.get_child("grouping");
		if (node.empty()) {
			setGroupingParams(node.get_value<std::string>());
		}
		else {
			setGroupingParams(asVector<int>(node));
		}
	}

	checkSetParam(tree, "sputtering", sputteringYield);

	if (tree.count("netParam")) {
		// Set the flag to not use the HDF5 file
		useHDF5Flag = false;
		auto node = tree.get_child("netParam");
		if (node.empty()) {
			setNetworkParameters(node.get_value<std::string>());
		}
		else {
			setNetworkParameters(asVector<IdType>(node));
		}
	}

	if (tree.count("radiusSize")) {
		auto node = tree.get_child("radiusSize");
		if (node.empty()) {
			setRadiusMinSizes(node.get_value<std::string>());
		}
		else {
			setRadiusMinSizes(asVector<int>(node));
		}
	}

	if (tree.count("boundary")) {
		auto node = tree.get_child("boundary");
		if (node.empty()) {
			setBoundaries(node.get_value<std::string>());
		}
		else {
			setBoundaries(asVector<int>(node));
		}
	}

	checkSetParam(tree, "xBCType", xBC);

	checkSetParam(tree, "heatLossPortion", heatLossPortion);

	checkSetParam(tree, "burstingDepth", burstingDepth);

	checkSetParam(tree, "burstingFactor", burstingFactor);

	if (tree.count("rng")) {
		processRNGParam(tree.get<std::string>("rng"));
	}

	checkSetParam(tree, "density", density);

	if (tree.count("pulse")) {
		setPulseParams(tree.get<std::string>("pulse"));
	}

	checkSetParam(tree, "lattice", latticeParameter);

	checkSetParam(tree, "impurityRadius", impurityRadius);

	checkSetParam(tree, "biasFactor", biasFactor);

	checkSetParam(tree, "hydrogenFactor", hydrogenFactor);

	checkSetParam(tree, "xenonDiffusivity", xenonDiffusivity);

	checkSetParam(tree, "fissionYield", fissionYield);

	checkSetParam(tree, "heVRatio", heVRatio);

	checkSetParam(tree, "migrationThreshold", migrationThreshold);

	checkSetParam(tree, "fluxDepthProfileFilePath", fluxDepthProfileFilePath);

	checkSetParam(tree, "basalPortion", basalPortion);

	checkSetParam(tree, "transitionSize", transitionSize);

	checkSetParam(tree, "cascadeDose", cascadeDose);

	checkSetParam(tree, "cascadeEfficiency", cascadeEfficiency);

	if (tree.count("petscArgs")) {
		for (auto&& elem : tree.get_child("petscArgs")) {
			appendPetscArg(elem.second.data());
		}
	}
}
} // namespace options
} // namespace xolotl
