#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <thread>

#include <boost/dll.hpp>
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <xolotl/interface/NetworkHandlerClassGenerator.h>
#include <xolotl/util/Filesystem.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/StreamUtils.h>
#include <xolotl/util/Tokenizer.h>

namespace xolotl
{
namespace interface
{
inline constexpr double nan = std::numeric_limits<double>::quiet_NaN();
inline constexpr AmountType maxAmount = std::numeric_limits<AmountType>::max();

std::ofstream
openFile(const fs::path& filePath)
{
	auto ofs = std::ofstream(filePath);
	if (!ofs) {
		XOLOTL_ERROR(
			std::runtime_error, "Unable to open file: " + filePath.string());
	}
	return ofs;
}

std::string
toLower(std::string s)
{
	std::transform(
		begin(s), end(s), begin(s), [](auto&& c) { return std::tolower(c); });
	return s;
}

std::string
removeSpaces(std::string s)
{
	s.erase(std::remove_if(
				begin(s), end(s), [](auto c) { return std::isspace(c); }),
		end(s));
	return s;
}

template <typename T, typename F>
void
printCommaList(std::ostream& os, const std::vector<T>& vec, F f)
{
	std::ostringstream oss;
	for (auto&& elem : vec) {
		oss << f(elem) << ", ";
	}
	auto str = oss.str();
	auto view = std::string_view(str);
	view.remove_suffix(2);
	os << view;
}

std::string
findReplace(std::string str, const std::string& from, const std::string& to)
{
	if (from.empty()) {
		return str;
	}
	std::size_t start_pos = 0;
	while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return str;
}

NetworkHandlerClassGenerator::NetworkHandlerClassGenerator(
	const fs::path& rnFile, const std::string& material,
	const std::string& userMaterial) :
	_rnFile(rnFile),
	_material(material),
	_userMaterial(userMaterial),
	_genDir(fs::current_path() / ".xgrn"),
	_buildDir(_genDir / "build"),
	_execFile(boost::dll::program_location().parent_path().string())
{
	readNetworkFile();
	_baseName = removeSpaces(_networkData.label) + "_u_";
	fs::create_directory(_genDir);
	_speciesList = _baseName + "SpeciesList";
	_prodReaction = _baseName + "ProductionReaction";
	_dissReaction = _baseName + "DissociationReaction";
	_sinkReaction = _baseName + "SinkReaction";
	_clusterGenerator = _baseName + "ClusterGenerator";
	_reactionGenerator = _baseName + "ReactionGenerator";
	_reactionNetwork = _baseName + "ReactionNetwork";
	_networkHandler = _baseName + "NetworkHandler";
	_materialHandler = _baseName + "MaterialHandler";

	_networkLibFile = _buildDir / ("lib" + _networkHandler + ".so");
}

void
NetworkHandlerClassGenerator::readNetworkFile()
{
	std::ifstream ifs(_rnFile);
	if (!ifs) {
		XOLOTL_ERROR(std::runtime_error,
			"Unable to open reaction network file: " + _rnFile.string());
	}

	auto ss = util::stripComments(ifs);
	boost::property_tree::iptree tree;
	boost::property_tree::read_json(ss, tree);

	// network
	auto nwNode = tree.get_child("network");
	_networkData.label = nwNode.get<std::string>("label");
	_networkData.desc = nwNode.get("description", "");
	if (!nwNode.count("lattice_parameter")) {
		XOLOTL_ERROR(std::invalid_argument,
			"\"lattice_parameter\" required in \"network\" group");
	}
	_networkData.latticeParameter = nwNode.get<double>("lattice_parameter");
	if (!nwNode.count("atomic_volume")) {
		XOLOTL_ERROR(std::invalid_argument,
			"\"atomic_volume\" required in \"network\" group");
	}
	_networkData.atomicVolume = nwNode.get<double>("atomic_volume");

	auto mapSpeciesType = [](const std::string& type) {
		if (type == "impurity") {
			return SpeciesData::Type::impurity;
		}
		else if (type == "vacancy") {
			return SpeciesData::Type::vacancy;
		}
		else if (type == "interstitial") {
			return SpeciesData::Type::interstitial;
		}
		else {
			XOLOTL_ERROR(std::invalid_argument,
				"unsupported species type: \"" + type + "\"");
		}
	};

	// species
	auto spNode = tree.get_child("species");
	for (auto&& [_, elNode] : spNode) {
		auto& spData = _speciesData.emplace_back();
		spData.name = elNode.get<std::string>("name");
		spData.label = elNode.get<std::string>("label");
		spData.type = mapSpeciesType(toLower(elNode.get<std::string>("type")));
		auto bIter = begin(elNode.get_child("bounds"));
		spData.bounds[0] = bIter->second.get_value<AmountType>();
		spData.bounds[1] = std::next(bIter)->second.get_value<AmountType>();
	}
	for (auto&& sp : _speciesData) {
		if (!_speciesLabelMap.try_emplace(sp.label, &sp).second) {
			XOLOTL_ERROR(std::invalid_argument,
				"Duplicate species label \"" + sp.label + "\"");
		}
	}

	auto checkSpeciesLabel = [*this](const std::string& label,
								 const std::string& msgPrefix) {
		auto it = std::find_if(begin(_speciesData), end(_speciesData),
			[&label](auto&& sp) { return sp.label == label; });
		if (it == end(_speciesData)) {
			XOLOTL_ERROR(std::invalid_argument,
				msgPrefix + " species (" + label +
					") not found in species list");
		}
	};

	// clusters
	auto clNode = tree.get_child("clusters");
	for (auto&& [_, elNode] : clNode) {
		auto& clData = _clusterData.emplace_back();
		clData.type = elNode.get<std::string>("type");
		checkSpeciesLabel(clData.type, "cluster");
		clData.size = {maxAmount, maxAmount};
		if (elNode.count("size")) {
			auto sizeNode = elNode.get_child("size");
			if (sizeNode.empty()) {
				clData.size[0] = sizeNode.get_value<AmountType>();
			}
			else {
				if (sizeNode.size() > 2) {
					XOLOTL_ERROR(std::invalid_argument,
						"cluster size must specify at most two values");
				}
				clData.size[0] =
					sizeNode.front().second.get_value<AmountType>();
				clData.size[1] = sizeNode.back().second.get_value<AmountType>();
			}
		}
		if (!elNode.count("radius")) {
			XOLOTL_ERROR(std::invalid_argument,
				"every cluster record must specify \"radius\"");
		}
		clData.radiusExpr = elNode.get<std::string>("radius");
		clData.migrationEnergy = elNode.get("migration_energy", nan);
		clData.diffusionFactor = elNode.get("diffusion_factor", nan);
	}
	// group clusters by species
	_clusterGroups.assign(_speciesData.size(), {});
	for (std::size_t i = 0; i < _speciesData.size(); ++i) {
		const ClusterData* general = nullptr;
		for (const auto& cl : _clusterData) {
			if (cl.type == _speciesData[i].label) {
				if (cl.size.size() == 0) {
					general = &cl;
				}
				else {
					_clusterGroups[i].push_back(&cl);
				}
			}
		}
		if (general != nullptr) {
			_clusterGroups[i].push_back(general);
		}
	}

	auto parseReactionExpr = [checkSpeciesLabel](const std::string& expr) {
		auto lr = util::Tokenizer<>{expr, "->"}();
		auto reactants = util::Tokenizer<>{lr[0], "+"}();
		auto products = util::Tokenizer<>{lr[1], "+"}();
		for (auto& r : reactants) {
			r = removeSpaces(r);
			checkSpeciesLabel(r, "reactant");
		}
		for (auto& p : products) {
			p = removeSpaces(p);
			if (p == "0") {
				if (products.size() > 1) {
					XOLOTL_ERROR(std::invalid_argument,
						"\"0\" can be used as reaction product only by itself");
				}
				products.clear();
				break;
			}
			checkSpeciesLabel(p, "product");
		}
		return std::make_tuple(reactants, products);
	};

	// reactions
	auto rnNode = tree.get_child("reactions");
	for (auto&& [_, elNode] : rnNode) {
		auto& rnData = _reactionData.emplace_back();
		if (elNode.count("expr")) {
			std::tie(rnData.reactantLabels, rnData.productLabels) =
				parseReactionExpr(elNode.get<std::string>("expr"));
		}
		else {
			if (!(elNode.count("reactants") && elNode.count("products"))) {
				XOLOTL_ERROR(std::runtime_error,
					"Reaction must be expressed either as \"expr\" or as "
					"\"reactants\" and \"products\"");
			}
			auto reactantNode = elNode.get_child("reactants");
			if (reactantNode.empty()) {
				XOLOTL_ERROR(std::runtime_error, "Reactant list empty");
			}
			for (auto&& reac : reactantNode) {
				auto s = reac.second.get_value<std::string>();
				checkSpeciesLabel(s, "reactant");
				rnData.reactantLabels.push_back(s);
			}

			auto productNode = elNode.get_child("products");
			if (productNode.empty()) {
				XOLOTL_ERROR(std::runtime_error, "Product list empty");
			}
			for (auto&& prod : productNode) {
				auto s = prod.second.get_value<std::string>();
				checkSpeciesLabel(s, "product");
				rnData.productLabels.push_back(s);
			}
		}
		for (auto&& r : rnData.reactantLabels) {
			rnData.reactants.push_back(_speciesLabelMap.at(r));
		}
		for (auto&& p : rnData.productLabels) {
			rnData.products.push_back(_speciesLabelMap.at(p));
		}
		rnData.type = elNode.get<std::string>("type");
		if (elNode.count("binding")) {
			rnData.bindingExpr = elNode.get<std::string>("binding");
		}
	}
	// group reactions by species set
	for (const auto& rn : _reactionData) {
		if (rn.reactants.size() > 1) {
			std::string key{};
			for (const auto& r : rn.reactantLabels) {
				key += r;
			}
			for (const auto& p : rn.productLabels) {
				key += p;
			}
			auto& grp = _productionReactionGroups[key];
			grp.general = &rn;
		}
		else {
			std::string key = rn.reactantLabels[0];
			for (const auto& prod : rn.productLabels) {
				key += prod;
			}
			auto& grp = _dissociationReactionGroups[key];
			grp.general = &rn;
		}
	}
}

void
NetworkHandlerClassGenerator::generateClusterGeneratorHeader()
{
	auto filePath = _genDir / (_clusterGenerator + ".h");
	auto ofs = openFile(filePath);
	ofs << "#pragma once\n";
	ofs << "#include <plsm/refine/Detector.h>\n";
	ofs << "#include <" << _baseName << "Traits.h>\n";

	ofs << "namespace xolotl::core::network {\n";

	auto detector = "plsm::refine::Detector<" + _clusterGenerator + ">";
	ofs << "class " << _clusterGenerator << " \n"
		<< "  : public " << detector << " {\n"
		<< "public:\n"
		<< "  using Species = " << _speciesList << ";\n"
		<< "  using Superclass = " << detector << ";\n"
		<< "  using NetworkType = " << _reactionNetwork << ";\n"
		<< "\n"
		   "  template <typename PlsmContext>\n"
		   "  using Cluster = typename NetworkType::Cluster<PlsmContext>;\n"
		   "\n"
		   "  using Region = typename NetworkType::Region;\n"
		   "  using Composition = typename NetworkType::Composition;\n"
		   "  using AmountType = typename NetworkType::AmountType;\n"
		   "  using BoolArray = plsm::refine::BoolVec<Region>;\n"
		   "\n"
		<< "  " << _clusterGenerator << "(const options::IOptions& options);\n"
		<< "\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  bool\n"
		   "  refine(const Region& region, BoolArray& result) const;\n"
		   "\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  bool\n"
		   "  select(const Region& region) const;\n"
		   "\n"
		   "  template <typename PlsmContext>\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  double\n"
		   "  getFormationEnergy(\n"
		   "    const Cluster<PlsmContext>& cluster) const noexcept;\n"
		   "\n"
		   "  template <typename PlsmContext>\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  double\n"
		   "  getMigrationEnergy(\n"
		   "    const Cluster<PlsmContext>& cluster) const noexcept;\n"
		   "\n"
		   "  template <typename PlsmContext>\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  double\n"
		   "  getDiffusionFactor(const Cluster<PlsmContext>& cluster,\n"
		   "    double latticeParameter) const noexcept;\n"
		   "\n"
		   "  template <typename PlsmContext>\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  double\n"
		   "  getReactionRadius(const Cluster<PlsmContext>& cluster,\n"
		   "    double latticeParameter, double interstitialBias,\n"
		   "    double impurityRadius) const noexcept;\n"
		   "\n"
		   "private:\n";
	ofs << "  AmountType _maxV{0};\n";
	ofs << "};\n";

	// close namespace
	ofs << "}\n";
}

void
NetworkHandlerClassGenerator::generateClusterGeneratorImpl()
{
	auto filePath = _genDir / (_clusterGenerator + ".tpp");
	auto ofs = openFile(filePath);
	ofs << "#pragma once\n"
		   "#include <xolotl/core/Constants.h>\n"
		   "#include <xolotl/util/MathUtils.h>\n"
		   "\n"
		   "namespace xolotl::core::network {\n"
		   "KOKKOS_INLINE_FUNCTION\n"
		   "bool\n"
		<< _clusterGenerator << "::refine(\n"
		<< "  const Region& region, BoolArray& result) const\n"
		   "{\n"
		   "  // no grouping\n"
		   "  result[0] = true;\n"
		   "  result[1] = true;\n"
		   "  return true;\n"
		   "}\n"
		   "\n"
		   "KOKKOS_INLINE_FUNCTION\n"
		   "bool\n"
		<< _clusterGenerator << "::select(const Region& region) const\n"
		<< "{\n"
        << "  int nAxis = 0;\n";
    for (auto&& spec : _speciesData) {
        ofs << "  nAxis += (region[Species::" << spec.label << "].begin() > 0);\n";
    }
    ofs << "  if (nAxis != 1) {\n"
           "    return false;\n"
           "  }\n"
		   "  return true;\n"
		   "}\n"
		   "\n"
		   "template <typename PlsmContext>\n"
		   "KOKKOS_INLINE_FUNCTION\n"
		   "double\n"
		<< _clusterGenerator << "::getFormationEnergy(\n"
		<< "  const Cluster<PlsmContext>& cluster) const noexcept\n"
		   "{\n"
		   "  return 0.0;\n" // TODO: need custom when using for binding energy
		   "}\n"
		   "\n";

	ofs << "template <typename PlsmContext>\n"
		   "KOKKOS_INLINE_FUNCTION\n"
		   "double\n"
		<< _clusterGenerator << "::getMigrationEnergy(\n"
		<< "  const Cluster<PlsmContext>& cluster) const noexcept\n"
		   "{\n"
		   "  const auto& reg = cluster.getRegion();\n"
		   "  if (reg.isSimplex()) {\n";
	ofs << "    Composition comp(reg.getOrigin());\n";
	for (auto&& clGroup : _clusterGroups) {
		if (clGroup.empty()) {
			continue;
		}
		auto species = "Species::" + clGroup[0]->type;
		ofs << "    if (comp.isOnAxis(" << species << ")) {\n"
			<< "      auto amt = comp[" << species << "];\n";
		for (auto cl : clGroup) {
			if (std::isnan(cl->migrationEnergy)) {
				continue;
			}
			if (cl->size[0] < maxAmount) {
				if (cl->size[1] < maxAmount) {
					ofs << "      if (" << cl->size[0] << " <= amt &&\n"
						<< "          amt <= " << cl->size[1] << ") {\n"
						<< "        return " << cl->migrationEnergy << ";\n"
						<< "      }\n";
				}
				else {
					ofs << "      if (amt == " << cl->size[0] << ") {\n"
						<< "        return " << cl->migrationEnergy << ";\n"
						<< "      }\n";
				}
			}
			else {
				ofs << "      return " << cl->migrationEnergy << ";\n";
			}
		}
		ofs << "    }\n";
	}
	ofs << "  }\n";
	ofs << "  return util::infinity<double>;\n"
		   "}\n"
		   "\n";

	ofs << "template <typename PlsmContext>\n"
		   "KOKKOS_INLINE_FUNCTION\n"
		   "double\n"
		<< _clusterGenerator << "::getDiffusionFactor(\n"
		<< "  const Cluster<PlsmContext>& cluster,\n"
		   "  double latticeParameter) const noexcept\n"
		   "{\n"
		   "  const auto& reg = cluster.getRegion();\n"
		   "  if (reg.isSimplex()) {\n";
	ofs << "    Composition comp(reg.getOrigin());\n";
	for (auto&& clGroup : _clusterGroups) {
		if (clGroup.empty()) {
			continue;
		}
		auto species = "Species::" + clGroup[0]->type;
		ofs << "    if (comp.isOnAxis(" << species << ")) {\n"
			<< "      auto amt = comp[" << species << "];\n";
		for (auto cl : clGroup) {
			if (std::isnan(cl->diffusionFactor)) {
				continue;
			}
			if (cl->size[0] < maxAmount) {
				if (cl->size[1] < maxAmount) {
					ofs << "      if (" << cl->size[0] << " <= amt &&\n"
						<< "          amt <= " << cl->size[1] << ") {\n"
						<< "        return " << cl->diffusionFactor << ";\n"
						<< "      }\n";
				}
				else {
					ofs << "      if (amt == " << cl->size[0] << ") {\n"
						<< "        return " << cl->diffusionFactor << ";\n"
						<< "      }\n";
				}
			}
			else {
				ofs << "      return " << cl->diffusionFactor << ";\n";
			}
		}
		ofs << "    }\n";
	}
	ofs << "  }\n"
		   "  return 0.0;\n"
		   "}\n"
		   "\n";

	const std::unordered_map<std::string, std::string> radiusKWMap{
		{"size", "amt"},
		{"atomic_volume", std::to_string(_networkData.atomicVolume)},
		{"pi", "::xolotl::core::pi"}, {"lattice_parameter", "latticeParameter"},
		{"impurity_radius", "impurityRadius"},
		{"interstitial_bias", "interstitialBias"}};
	auto parseRadiusExpr = [&radiusKWMap](std::string expr) {
		for (auto [from, to] : radiusKWMap) {
			expr = findReplace(expr, from, to);
		}
		return expr;
	};

	ofs << "template <typename PlsmContext>\n"
		   "KOKKOS_INLINE_FUNCTION\n"
		   "double\n"
		<< _clusterGenerator << "::getReactionRadius(\n"
		<< "  const Cluster<PlsmContext>& cluster,\n"
		   "  double latticeParameter, double interstitialBias,\n"
		   "  double impurityRadius) const noexcept\n"
		   "{\n"
		   "  const auto& reg = cluster.getRegion();\n"
		   "  if (reg.isSimplex()) {\n";
	ofs << "    Composition comp(reg.getOrigin());\n";
	for (auto&& clGroup : _clusterGroups) {
		if (clGroup.empty()) {
			continue;
		}
		auto species = "Species::" + clGroup[0]->type;
		ofs << "    if (comp.isOnAxis(" << species << ")) {\n"
			<< "      auto amt = comp[" << species << "];\n";
		for (auto cl : clGroup) {
			const auto& expr = cl->radiusExpr;
			if (cl->size[0] < maxAmount) {
				if (cl->size[1] < maxAmount) {
					ofs << "      if (" << cl->size[0] << " <= amt &&\n"
						<< "          amt <= " << cl->size[1] << ") {\n"
						<< "        // parsed from: " << expr << "\n"
						<< "        return " << parseRadiusExpr(expr) << ";\n"
						<< "      }\n";
				}
				else {
					ofs << "      if (amt == " << cl->size[0] << ") {\n"
						<< "        // parsed from: " << expr << "\n"
						<< "        return " << parseRadiusExpr(expr) << ";\n"
						<< "      }\n";
				}
			}
			else {
				ofs << "      // parsed from: " << expr << "\n"
					<< "      return " << parseRadiusExpr(expr) << ";\n";
			}
		}
		ofs << "    }\n";
	}
	ofs << "  }\n"
		   "  return 0.0;\n"
		   "}\n"
		   "}\n";

	filePath = _genDir / (_clusterGenerator + ".cpp");
	ofs = openFile(filePath);
	ofs << "#include <" << _reactionNetwork << ".h>\n"
		<< "namespace xolotl::core::network {\n"
		<< _clusterGenerator << "::" << _clusterGenerator << "(\n"
		<< "  const options::IOptions& options) :\n"
		   "  _maxV(options.getMaxV())\n"
		   "{\n"
		   "}\n"
		   "}\n";
}

void
NetworkHandlerClassGenerator::generateClusterGenerator()
{
	generateClusterGeneratorHeader();
	generateClusterGeneratorImpl();
}

void
NetworkHandlerClassGenerator::generateReactionHeader()
{
	auto filePath = _genDir / (_baseName + "Reaction.h");
	auto ofs = openFile(filePath);
	ofs << "#pragma once\n";
	ofs << "#include <" << _baseName << "Traits.h>\n";
	ofs << "#include <xolotl/core/network/SinkReaction.h>\n";

	ofs << "namespace xolotl::core::network {\n";

	ofs << "class " << _reactionNetwork << ";\n\n";

	ofs << "class " << _prodReaction << " :\n"
		<< "  public ProductionReaction<" << _reactionNetwork << ", "
		<< _prodReaction << ">\n";
	ofs << "{\n" << "public:\n";
	ofs << "  using Superclass = ProductionReaction<" << _reactionNetwork
		<< ", " << _prodReaction << ">;\n";
	ofs << "  using Superclass::Superclass;\n"
		   "\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  double\n"
		   "  getRateForProduction(IndexType gridIndex);\n";
	ofs << "};\n\n";

	ofs << "class " << _dissReaction << " :\n"
		<< "  public DissociationReaction<" << _reactionNetwork << ", "
		<< _dissReaction << ">\n";
	ofs << "{\n" << "public:\n";
	ofs << "  using Superclass = DissociationReaction<" << _reactionNetwork
		<< ", " << _dissReaction << ">;\n";
	ofs << "  using Superclass::Superclass;\n"
		   "\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  double\n"
		   "  getRateForProduction(IndexType gridIndex);\n"
		   "\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  double\n"
		   "  computeBindingEnergy(double time = 0.0);\n";
	ofs << "};\n\n";

	ofs << "class " << _sinkReaction << " :\n"
		<< "  public SinkReaction<" << _reactionNetwork << ", " << _sinkReaction
		<< ">\n";
	ofs << "{\n" << "public:\n";
	ofs << "  using Superclass = SinkReaction<" << _reactionNetwork << ", "
		<< _sinkReaction << ">;\n";
	ofs << "  using Superclass::Superclass;\n"
		   "\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  double\n"
		   "  getSinkBias();\n"
		   "\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  double\n"
		   "  getSinkStrength();\n";
	ofs << "};\n";

	ofs << "}\n";
}

void
NetworkHandlerClassGenerator::generateReactionImpl()
{
	auto filePath = _genDir / (_baseName + "Reaction.tpp");
	auto ofs = openFile(filePath);
	ofs << "#pragma once\n";

	ofs << "#include <xolotl/core/network/impl/SinkReaction.tpp>\n";
	ofs << "#include <xolotl/util/MathUtils.h>\n";
	ofs << "namespace xolotl::core::network {\n";

	ofs << "namespace " + _baseName + " {\n";
	ofs << "  template <typename TRegion>\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  double\n"
		   "  getRate(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg, \n"
		   "      const double r0, const double r1, const double dc0, \n"
		   "      const double dc1) {\n"
		   "    constexpr double pi = ::xolotl::core::pi;\n"
		   "    double kPlus = 4.0 * pi * (r0 + r1) * (dc0 + dc1);\n"
		   "    return kPlus;\n"
		   "  }\n";
	ofs << "}\n";

	ofs << "KOKKOS_INLINE_FUNCTION\n"
		   "double\n"
		<< _prodReaction << "::getRateForProduction(IndexType gridIndex) {\n"
		<< "  auto cl0 = this->_clusterData->getCluster(_reactants[0]);\n"
		   "  auto cl1 = this->_clusterData->getCluster(_reactants[1]);\n"
		   "  double r0 = cl0.getReactionRadius();\n"
		   "  double r1 = cl1.getReactionRadius();\n"
		   "  double dc0 = cl0.getDiffusionCoefficient(gridIndex);\n"
		   "  double dc1 = cl1.getDiffusionCoefficient(gridIndex);\n"
		<< "  return " << _baseName << "::getRate(\n"
		<< "    cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);\n"
		   "}\n\n";

	ofs << "KOKKOS_INLINE_FUNCTION\n"
		   "double\n"
		<< _dissReaction << "::getRateForProduction(IndexType gridIndex) {\n"
		<< "  auto cl0 = this->_clusterData->getCluster(_products[0]);\n"
		   "  auto cl1 = this->_clusterData->getCluster(_products[1]);\n"
		   "  double r0 = cl0.getReactionRadius();\n"
		   "  double r1 = cl1.getReactionRadius();\n"
		   "  double dc0 = cl0.getDiffusionCoefficient(gridIndex);\n"
		   "  double dc1 = cl1.getDiffusionCoefficient(gridIndex);\n"
		<< "  return " << _baseName << "::getRate(\n"
		<< "    cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);\n"
		   "}\n\n";

	const std::unordered_map<std::string, std::string> bindingKWMap{
		{"size", "amt"}};
	auto parseBindingExpr = [&bindingKWMap](std::string expr) {
		for (auto [from, to] : bindingKWMap) {
			expr = findReplace(expr, from, to);
		}
		return expr;
	};

	ofs << "KOKKOS_INLINE_FUNCTION\n"
		   "double\n"
		<< _dissReaction << "::computeBindingEnergy(double time) {\n"
		<< "  using Species = typename Superclass::Species;\n"
		   "  using Composition = typename Superclass::Composition;\n"
		   "  double be = 5.0;\n"
		   "  auto cl = this->_clusterData->getCluster(this->_reactant);\n"
		   "  auto prod1 =\n"
		   "    this->_clusterData->getCluster(this->_products[0]);\n"
		   "  auto prod2 =\n"
		   "    this->_clusterData->getCluster(this->_products[1]);\n"
		   "  auto clReg = cl.getRegion();\n"
		   "  auto prod1Reg = prod1.getRegion();\n"
		   "  auto prod2Reg = prod2.getRegion();\n"
		   "  if (clReg.isSimplex() && prod1Reg.isSimplex() &&\n"
		   "      prod2Reg.isSimplex()) {\n";
	for (auto&& [key, rnGroup] : _dissociationReactionGroups) {
		if (!rnGroup.general) {
			XOLOTL_ERROR(std::invalid_argument,
				"no general case provided for reaction " + key);
		}
		auto rSpec = rnGroup.general->reactantLabels[0];
		auto p1Spec = rnGroup.general->productLabels[0];
		auto p2Spec = rnGroup.general->productLabels[1];
		const auto& expr = rnGroup.general->bindingExpr;
		ofs << "    Composition rComp = clReg.getOrigin();\n"
			   "    Composition p1Comp = prod1Reg.getOrigin();\n"
			   "    Composition p2Comp = prod2Reg.getOrigin();\n"
			<< "    if (rComp.isOnAxis(Species::" << rSpec << ") &&\n"
			<< "        p1Comp.isOnAxis(Species::" << p1Spec << ") &&\n"
			<< "        p2Comp.isOnAxis(Species::" << p2Spec << ")) {\n"
			<< "      auto amt = rComp[Species::" << rSpec << "];\n";
		for (auto rn : rnGroup.cases) {
			// TODO: cases for specific sizes
		}
		ofs << "      // parsed from: " << expr << "\n"
			<< "      be = " << parseBindingExpr(expr) << ";\n"
			<< "    }\n";
	}
	ofs << "  }\n"
            " return util::max(-5.0, util::min(be, 5.0));\n"
		   // "  return util::clamp(be, -5.0, 5.0);\n"
		   "}\n\n";

	ofs << "KOKKOS_INLINE_FUNCTION\n"
		   "double\n"
		<< _sinkReaction << "::getSinkBias() {\n"
		<< "  return 1.0;\n" // TODO: need for general case
		   "}\n\n";

	ofs << "KOKKOS_INLINE_FUNCTION\n"
		   "double\n"
		<< _sinkReaction << "::getSinkStrength() {\n"
		<< "  return 0.0;\n" // TODO: need for general case
		   "}\n";

	// close namespace
	ofs << "}\n";
}

void
NetworkHandlerClassGenerator::generateReactions()
{
	generateReactionHeader();
	generateReactionImpl();
}

void
NetworkHandlerClassGenerator::generateTraits()
{
	auto filePath = _genDir / (_baseName + "Traits.h");
	auto ofs = openFile(filePath);
	ofs << "#pragma once\n"
		   "#include <tuple>\n"
		   "#include <xolotl/core/network/ReactionNetworkTraits.h>\n";

	ofs << "namespace xolotl::core::network {\n";

	ofs << "class " << _reactionNetwork << ";\n";
	ofs << "class " << _prodReaction << ";\n";
	ofs << "class " << _dissReaction << ";\n";
	ofs << "class " << _sinkReaction << ";\n";
	ofs << "class " << _clusterGenerator << ";\n\n";

	ofs << "enum class " << _speciesList << " { ";
	printCommaList(ofs, _speciesData, [](auto&& s) { return s.label; });
	ofs << " };\n\n";

	ofs << "inline const std::string&"
		<< "toLabelString(" << _speciesList << " species) {\n";
	ofs << "  static const std::string labelArray[] = { ";
	printCommaList(
		ofs, _speciesData, [](auto&& s) { return std::quoted(s.label); });
	ofs << " };\n";
	ofs << "  return labelArray[static_cast<int>(species)];\n";
	ofs << "}\n\n";

	ofs << "inline const std::string&"
		<< "toNameString(" << _speciesList << " species) {\n";
	ofs << "  static const std::string nameArray[] = { ";
	printCommaList(
		ofs, _speciesData, [](auto&& s) { return std::quoted(s.name); });
	ofs << " };\n";
	ofs << "  return nameArray[static_cast<int>(species)];\n";
	ofs << "}\n\n";

	auto nSpecies = _speciesData.size();
	ofs << "template <>\n";
	ofs << "struct NumberOfSpecies<" << _speciesList << "> : \n"
		<< "  std::integral_constant<std::size_t, " << nSpecies << "> {};\n\n";

	auto nInterstitial = std::count_if(begin(_speciesData), end(_speciesData),
		[](auto&& sp) { return sp.type == SpeciesData::Type::interstitial; });
	ofs << "template <>\n";
	ofs << "struct NumberOfInterstitialSpecies<" << _speciesList << "> : \n"
		<< "  std::integral_constant<std::size_t, " << nInterstitial
		<< "> {};\n\n";

	auto nVacancy = std::count_if(begin(_speciesData), end(_speciesData),
		[](auto&& sp) { return sp.type == SpeciesData::Type::vacancy; });
	ofs << "template <>\n";
	ofs << "struct NumberOfVacancySpecies<" << _speciesList << "> : \n"
		<< "  std::integral_constant<std::size_t, " << nVacancy << "> {};\n\n";

	// TODO: SpeciesForGrouping (default is fine for now)

	ofs << "template <>\n";
	ofs << "struct ReactionNetworkTraits<" << _baseName
		<< "ReactionNetwork> {\n";
	ofs << "  using Species = " << _speciesList << ";\n";
	ofs << "  static constexpr std::size_t numSpecies = " << nSpecies << ";\n";
	ofs << "  using ProductionReactionType = " << _prodReaction << ";\n";
	ofs << "  using DissociationReactionType = " << _dissReaction << ";\n";
	ofs << "  using SinkReactionType = " << _sinkReaction << ";\n";

	ofs << "  using ReactionTypeList = std::tuple<ProductionReactionType,\n"
		   "    DissociationReactionType, SinkReactionType>;\n";

	ofs << "  using ClusterGenerator = " << _clusterGenerator << ";\n";
	ofs << "};\n";

	// close namespace
	ofs << "}\n";
}

void
NetworkHandlerClassGenerator::generateNetworkHeader()
{
	auto filePath = _genDir / (_reactionNetwork + ".h");
	auto ofs = openFile(filePath);
	ofs << "#pragma once\n";
	ofs << "#include <" << _baseName << "Reaction.h>\n";
	ofs << "#include <" << _baseName << "Traits.h>\n";
	ofs << "#include <xolotl/core/network/ReactionNetwork.h>\n";

	ofs << "namespace xolotl::core::network {\n";

	ofs << "namespace detail\n"
		<< "{\n"
		<< "class " << _reactionGenerator << ";\n"
		<< "}\n";

	ofs << "class " << _reactionNetwork << "\n"
		<< "  : public ReactionNetwork<" << _reactionNetwork << ">\n"
		<< "{\n"
		<< "  friend class ReactionNetwork<" << _reactionNetwork << ">;\n"
		<< "public:\n"
		<< "  using Superclass = ReactionNetwork<" << _reactionNetwork << ">;\n"
		<< "  using Subpaving = typename Superclass::Subpaving;\n"
		   "  using Composition = typename Superclass::Composition;\n"
		   "  using Species = typename Superclass::Species;\n"
		   "  using AmountType = typename Superclass::AmountType;\n"
		   "  using IndexType = typename Superclass::IndexType;\n"
		   "  using ConcentrationsView =\n"
		   "    typename Superclass::ConcentrationsView;\n"
		   "  using FluxesView = typename Superclass::FluxesView;\n"
		   "\n"
		   "  using Superclass::Superclass;\n"
		   "\n"
		   "  IndexType\n"
		   "  checkLargestClusterId();\n"
		   "\n"
		   "  std::string\n"
		   "  getMonitorOutputFileName() const override\n"
		   "  {\n"
		   "    return \"monitor.dat\";\n"
		   "  }\n"
		   "\n"
		   "  std::string\n"
		   "  getMonitorDataHeaderString() const override;\n"
		   "\n"
		   "  void\n"
		   "  addMonitorDataValues(Kokkos::View<const double*> conc,\n"
		   "    double fac, std::vector<double>& totalVals) override;\n"
		   "\n"
		   "  std::size_t\n"
		   "  getMonitorDataLineSize() const override;\n"
		   "\n"
		   "  void\n"
		   "  writeMonitorDataLine(const std::vector<double>& localData,\n"
		   "    double time) override;\n"
		   "\n"
		   "private:\n"
		   "  double\n"
		   "  checkLatticeParameter(double latticeParameter);\n"
		   "\n"
		   "  double\n"
		   "  computeAtomicVolume(double latticeParameter);\n"
		   "\n"
		   "  double\n"
		   "  checkImpurityRadius(double impurityRadius);\n"
		   "\n"
		<< "  detail::" << _reactionGenerator << "\n"
		<< "  getReactionGenerator() const noexcept;\n"
		   "\n"
		   "  void\n"
		   "  readClusters(const std::string&) { }\n"
		   "\n"
		   "  void\n"
		   "  readReactions(double, const std::string&) { }\n"
		   "};\n\n";

	ofs << "namespace detail {\n"
		<< "class " << _reactionGenerator << "\n"
		<< "  : public ReactionGenerator<" << _reactionNetwork << ", "
		<< _reactionGenerator << ">\n"
		<< "{\n"
		<< "  friend class ReactionGeneratorBase<" << _reactionNetwork << ", "
		<< _reactionGenerator << ">;\n"
		<< "\n"
		<< "public:\n"
		<< "  using NetworkType = " << _reactionNetwork << ";\n"
		<< "  using Subpaving = typename NetworkType::Subpaving;\n"
		   "  using IndexType = typename NetworkType::IndexType;\n"
		   "\n"
		   "  using Superclass = \n"
		<< "    ReactionGenerator<" << _reactionNetwork << ", "
		<< _reactionGenerator << ">;\n"
		<< "\n"
		   "  using Superclass::Superclass;\n"
		   "\n"
		   "  template <typename TTag>\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  void\n"
		   "  operator()(IndexType i, IndexType j, TTag tag) const;\n"
		   "\n"
		   "  template <typename TTag>\n"
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  void\n"
		   "  addSinks(IndexType i, TTag tag) const;\n"
		   "\n"
		   "private:\n"
		   "  ReactionCollection<NetworkType>\n"
		   "  getReactionCollection() const;\n"
		   "};\n"
		   "} // namespace detail\n";

	// close namespace
	ofs << "}\n";

	ofs << "#include <" << _clusterGenerator << ".h>\n"
		<< "#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)\n"
		<< "#include <" << _reactionNetwork << ".tpp>\n"
		<< "#endif\n";
}

void
NetworkHandlerClassGenerator::generateNetworkImpl()
{
	auto filePath = _genDir / (_reactionNetwork + ".tpp");
	auto ofs = openFile(filePath);
	ofs << "#pragma once\n"
		   "#include <xolotl/core/Constants.h>\n"
		   "#include "
		   "<xolotl/core/network/detail/impl/SinkReactionGenerator.tpp>\n"
		<< "#include <" << _clusterGenerator << ".tpp>\n"
		<< "#include <" << _baseName << "Reaction.tpp>\n"
		<< "#include <xolotl/core/network/impl/ReactionNetwork.tpp>\n"
		   "\n"
		   "namespace xolotl::core::network\n"
		   "{\n"
		   "namespace detail\n"
		   "{\n"
		   "template <typename TTag>\n"
		   "KOKKOS_INLINE_FUNCTION\n"
		   "void\n"
		<< _reactionGenerator << "::operator()(\n"
		<< "  IndexType i, IndexType j, TTag tag) const\n"
		   "{\n"
		   "  auto diffusionFactor = this->_clusterData.diffusionFactor;\n"
		   "  if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {\n"
		   "    return;\n"
		   "  }\n"
		   "\n"
		   "  using Species = typename NetworkType::Species;\n"
		   "  using Composition = typename NetworkType::Composition;\n"
		   "  using AmountType = typename NetworkType::AmountType;\n"
		   "\n"
		   "  constexpr auto species = NetworkType::getSpeciesRange();\n"
		   "  constexpr auto speciesNoI = NetworkType::getSpeciesRangeNoI();\n"
		   "\n"
		   "  const auto& cl1Reg = this->getCluster(i).getRegion();\n"
		   "  const auto& cl2Reg = this->getCluster(j).getRegion();\n"
		   "  Composition lo1 = cl1Reg.getOrigin();\n"
		   "  Composition lo2 = cl2Reg.getOrigin();\n"
		   "\n"
		   "  auto& subpaving = this->getSubpaving();\n"
		   "\n"
		   "  if (cl1Reg.isSimplex() && cl2Reg.isSimplex()) {\n";

	auto checkVacItrPair = [](const SpeciesData* r1, const SpeciesData* r2) {
		constexpr auto vac = SpeciesData::Type::vacancy;
		constexpr auto itr = SpeciesData::Type::interstitial;
		if (r1->type == vac && r2->type == itr) {
			return std::make_tuple(true, r1, r2);
		}
		else if (r2->type == vac && r1->type == itr) {
			return std::make_tuple(true, r2, r1);
		}
		const SpeciesData* nul{nullptr};
		return std::make_tuple(false, nul, nul);
	};

	for (auto&& [key, grp] : _productionReactionGroups) {
		auto r1 = grp.general->reactants[0];
		auto r2 = grp.general->reactants[1];
		const auto& products = grp.general->products;
		if (r1->label == r2->label) {
			auto spec = r1->label;
			auto specRef = "Species::" + spec;
			auto pSpRef = "Species::" + products[0]->label;
			ofs << "    if (lo1.isOnAxis(" << specRef << ") &&\n"
				<< "        lo2.isOnAxis(" << specRef << ")) {\n"
				<< "      auto size = lo1[" << specRef << "] +\n"
				<< "        lo2[" << specRef << "];\n"
				<< "      Composition comp = Composition::zero();\n"
				<< "      comp[" << pSpRef << "] = size;\n"
				<< "      auto pId = subpaving.findTileId(comp);\n"
				<< "      if (pId != subpaving.invalidIndex()) {\n"
				<< "        this->addProductionReaction(tag, {i, j, pId});\n"
				<< "      }\n"
				<< "    }\n";
		}
		else {
			auto [flag, v, i] = checkVacItrPair(r1, r2);
			if (flag) {
				auto vSpec = v->label;
				auto iSpec = i->label;
				auto vSpecRef = "Species::" + vSpec;
				auto iSpecRef = "Species::" + iSpec;
				ofs << "    if ((lo1.isOnAxis(" << vSpecRef << ") &&\n"
					<< "        lo2.isOnAxis(" << iSpecRef << ")) ||\n"
					<< "        (lo1.isOnAxis(" << iSpecRef << ") &&\n"
					<< "        lo2.isOnAxis(" << vSpecRef << "))) {\n"
					<< "      auto vSize = lo1.isOnAxis(" << vSpecRef << ") ?\n"
					<< "          lo1[" << vSpecRef << "] :\n"
					<< "          lo2[" << vSpecRef << "];\n"
					<< "      auto iSize = lo1.isOnAxis(" << iSpecRef << ") ?\n"
					<< "          lo1[" << iSpecRef << "] :\n"
					<< "          lo2[" << iSpecRef << "];\n"
					<< "      auto prodSize = vSize - iSize;\n";
				if (products.empty()) {
					ofs << "      if (prodSize == 0) {\n"
						   "        this->addProductionReaction(tag, {i, j});\n"
						   "      }\n";
				}
				else if (products.size() == 1 &&
					products[0]->type == SpeciesData::Type::vacancy) {
					ofs << "      if (prodSize > 0) {\n"
						<< "        Composition comp = Composition::zero();\n"
						<< "        comp[" << vSpecRef << "] = prodSize;\n"
						<< "        auto pId = subpaving.findTileId(comp);\n"
						   "        if (pId != subpaving.invalidIndex()) {\n"
						   "          this->addProductionReaction(\n"
						   "            tag, {i, j, pId});\n"
						   "        }\n"
						   "      }\n";
				}
				else if (products.size() == 1 &&
					products[0]->type == SpeciesData::Type::interstitial) {
					ofs << "      if (prodSize < 0) {\n"
						<< "        Composition comp = Composition::zero();\n"
						<< "        comp[" << iSpecRef << "] = -prodSize;\n"
						<< "        auto pId = subpaving.findTileId(comp);\n"
						   "        if (pId != subpaving.invalidIndex()) {\n"
						   "          this->addProductionReaction(\n"
						   "            tag, {i, j, pId});\n"
						   "        }\n"
						   "      }\n";
				}
				else {
					XOLOTL_ERROR(
						std::runtime_error, "unsupported reaction case");
				}
				ofs << "    }\n";
			}
			else {
				XOLOTL_ERROR(std::runtime_error, "unsupported reaction case");
			}
		}
	}

	for (auto&& [key, grp] : _dissociationReactionGroups) {
		auto r1 = grp.general->products[0];
		auto r2 = grp.general->products[1];
		const auto& reactants = grp.general->reactants;
		if (r1->label == r2->label) {
			auto spec = r1->label;
			auto specRef = "Species::" + spec;
			auto rSpRef = "Species::" + reactants[0]->label;
			ofs << "    if (lo1.isOnAxis(" << specRef << ") &&\n"
				<< "        lo2.isOnAxis(" << specRef << ")) {\n"
				<< "      auto size = lo1[" << specRef << "] +\n"
				<< "        lo2[" << specRef << "];\n"
				<< "      Composition comp = Composition::zero();\n"
				<< "      comp[" << rSpRef << "] = size;\n"
				<< "      auto pId = subpaving.findTileId(comp);\n"
				<< "      if (pId != subpaving.invalidIndex()) {\n"
				<< "        this->addDissociationReaction(tag, {pId, i, j});\n"
				<< "      }\n"
				<< "    }\n";
		}
		else {
			XOLOTL_ERROR(std::runtime_error, "unsupported reaction case");
		}
	}
	ofs << "  }\n"
		   "}\n"
		   "\n"
		   "template <typename TTag>\n"
		   "KOKKOS_INLINE_FUNCTION\n"
		   "void\n"
		<< _reactionGenerator << "::addSinks(IndexType i, TTag tag) const\n"
		<< "{\n"
		   "" // TODO: no sinks for now
		   "}\n"
		   "\n"
		<< "inline\n"
		<< "ReactionCollection<" << _reactionGenerator << "::NetworkType>\n"
		<< _reactionGenerator << "::getReactionCollection() const\n"
		<< "{\n"
		   "  ReactionCollection<NetworkType> ret(\n"
		   "    this->_clusterData.gridSize, this->_clusterData.numClusters,\n"
		   "    this->_enableReadRates, this->getProductionReactions(),\n"
		   "    this->getDissociationReactions(), this->getSinkReactions());\n"
		   "  return ret;\n"
		   "}\n";

	ofs << "}\n"; // namespace detail

	ofs << "inline detail::" << _reactionGenerator << "\n"
		<< _reactionNetwork << "::getReactionGenerator() const noexcept\n"
		<< "{\n"
		<< "  return detail::" << _reactionGenerator << "{*this};\n"
		<< "}\n";

	ofs << "}\n";

	filePath = _genDir / (_reactionNetwork + ".cpp");
	ofs = openFile(filePath);
	ofs << "#include <" << _reactionNetwork << ".h>\n"
		<< "#include <" << _reactionNetwork << ".tpp>\n"
		<< "#include <xolotl/util/MPIUtils.h>\n"
		<< "\n"
		<< "namespace xolotl::core::network\n"
		<< "{\n"
		<< "template\n"
		<< "ReactionNetwork<" << _reactionNetwork << ">::ReactionNetwork();\n"
		<< "\n"
		<< "template\n"
		<< "ReactionNetwork<" << _reactionNetwork << ">::ReactionNetwork(\n"
		<< "  const std::vector<AmountType>& maxSpeciesAmounts,\n"
		   "  const std::vector<SubdivisionRatio>& subdivisionRatios,\n"
		   "  IndexType gridSize,\n"
		   "  const options::IOptions& opts);\n"
		   "\n"
		<< "template\n"
		<< "ReactionNetwork<" << _reactionNetwork << ">::ReactionNetwork(\n"
		<< "  const std::vector<AmountType>& maxSpeciesAmounts,\n"
		   "  IndexType gridSize,\n"
		   "  const options::IOptions& opts);\n"
		   "\n"
		   "template\n"
		<< "ReactionNetwork<" << _reactionNetwork << ">::~ReactionNetwork();\n"
		<< "\n"
		<< "template double\n"
		<< "ReactionNetwork<" << _reactionNetwork << ">\n"
		<< "  ::getTotalConcentration(ConcentrationsView concentrations,\n"
		   "    Species type, AmountType minSize);\n"
		<< "\n"
		<< "template double\n"
		<< "ReactionNetwork<" << _reactionNetwork << ">\n"
		<< "  ::getTotalRadiusConcentration(\n"
		   "    ConcentrationsView concentrations,\n"
		   "    Species type, AmountType minSize);\n"
		<< "\n"
		<< "template double\n"
		<< "ReactionNetwork<" << _reactionNetwork << ">\n"
		<< "  ::getTotalAtomConcentration(ConcentrationsView concentrations,\n"
		   "    Species type, AmountType minSize);\n"
		<< "\n"
		<< "template double\n"
		<< "ReactionNetwork<" << _reactionNetwork << ">\n"
		<< "  ::getTotalTrappedAtomConcentration(\n"
		   "    ConcentrationsView concentrations,\n"
		   "    Species type, AmountType minSize);\n"
		<< "\n"
		<< "template double\n"
		<< "ReactionNetwork<" << _reactionNetwork << ">\n"
		<< "  ::getTotalVolumeFraction(ConcentrationsView concentrations,\n"
		   "    Species type, AmountType minSize);\n"
		   "\n"
		   "double\n"
		<< _reactionNetwork
		<< "::checkLatticeParameter(double latticeParameter)\n"
		<< "{\n"
		<< "  constexpr double lpConstant = " << _networkData.latticeParameter
		<< ";\n"
		<< "  if (latticeParameter != lpConstant) {\n"
		   "    if (latticeParameter > 0.0) {\n"
		   "      XOLOTL_LOG_WARN\n"
		   "        << \"Changing lattice parameter in custom network\";\n"
		   "      return latticeParameter;\n"
		   "    }\n"
		   "  }\n"
		   "  return lpConstant;\n"
		   "}\n"
		   "\n"
		   "double\n"
		<< _reactionNetwork
		<< "::computeAtomicVolume(double latticeParameter)\n"
		<< "{\n"
		<< "  return " << _networkData.atomicVolume << ";\n"
		<< "}\n"
		   "\n"
		   "double\n"
		<< _reactionNetwork << "::checkImpurityRadius(double impurityRadius)\n"
		<< "{\n"
		   "  return impurityRadius;\n" // TODO: ?
		   "}\n"
		   "\n"
		<< _reactionNetwork << "::IndexType\n"
		<< _reactionNetwork << "::checkLargestClusterId()\n"
		<< "{\n"
		   "  auto clData = _clusterData.d_view;\n"
		   "  using Reducer = Kokkos::MaxLoc<AmountType, IndexType>;\n"
		   "  Reducer::value_type maxLoc;\n"
		   "  Kokkos::parallel_reduce(\n"
		   "    \""
		<< _reactionNetwork << "::checkLargestClusterId\",\n"
		<< "    _numClusters,\n"
		   "    KOKKOS_LAMBDA(IndexType i, Reducer::value_type& update) {\n"
		   "      const Region& clReg = clData().getCluster(i).getRegion();\n"
		   "      Composition hi = clReg.getUpperLimitPoint();\n"
		   "      auto size = hi[Species::V];\n" // TODO: ?
		   "      if (size > update.val) {\n"
		   "        update.val = size;\n"
		   "        update.loc = i;\n"
		   "      }\n"
		   "    },\n"
		   "    Reducer(maxLoc));\n"
		   "  return maxLoc.loc;\n"
		   "}\n"
		   "\n";

	std::size_t nVars = 3;

	ofs << "std::string\n"
		<< _reactionNetwork << "::getMonitorDataHeaderString() const\n"
		<< "{\n"
		   "  std::stringstream header;\n"
		   "  auto numSpecies = getSpeciesListSize();\n"
		   "  header << \"#time \";\n"
		   "  for (auto id = SpeciesId(numSpecies); id; ++id) {\n"
		   "    auto speciesName = this->getSpeciesName(id);\n"
		   "    header << speciesName << \"_density \"\n"
		   "      << speciesName << \"_atom \"\n"
		   "      << speciesName << \"_diameter \"\n"
		   "      << speciesName << \"_partial_density \"\n"
		   "      << speciesName << \"_partial_atom \"\n"
		   "      << speciesName << \"_partial_diameter \";\n"
		   "  }\n"
		   "  return header.str();\n"
		   "}\n"
		   "\n";

	ofs << "void\n"
		<< _reactionNetwork << "::addMonitorDataValues(\n"
		<< "  Kokkos::View<const double*> conc, double fac,\n"
		   "  std::vector<double>& totalVals)\n"
		   "{\n"
		   "  auto numSpecies = getSpeciesListSize();\n"
		   "  const auto& minSizes = this->getMinRadiusSizes();\n"
		   "  for (auto id = SpeciesId(numSpecies); id; ++id) {\n"
		   "    using TQ = IReactionNetwork::TotalQuantity;\n"
		   "    using Q = TQ::Type;\n"
		<< "    using TQA = util::Array<TQ, " << 2 * nVars << ">;\n"
		<< "    auto ms = minSizes[id()];\n"
		   "    auto totals = this->getTotals(conc,\n"
		   "      TQA{TQ{Q::total, id, 1},\n"
		   "        TQ{Q::atom, id, 1},\n"
		   "        TQ{Q::radius, id, 1},\n"
		   "        TQ{Q::total, id, ms},\n"
		   "        TQ{Q::atom, id, ms},\n"
		   "        TQ{Q::radius, id, ms}});\n"
		<< "    totalVals[(" << 2 * nVars
		<< " * id()) + 0] += totals[0] * fac;\n"
		<< "    totalVals[(" << 2 * nVars
		<< " * id()) + 1] += totals[1] * fac;\n"
		<< "    totalVals[(" << 2 * nVars
		<< " * id()) + 2] += totals[2] * fac;\n"
		<< "    totalVals[(" << 2 * nVars
		<< " * id()) + 3] += totals[3] * fac;\n"
		<< "    totalVals[(" << 2 * nVars
		<< " * id()) + 4] += totals[4] * fac;\n"
		<< "    totalVals[(" << 2 * nVars
		<< " * id()) + 5] += totals[5] * fac;\n"
		<< "\n"
		<< "    totalVals[(" << 2 * nVars << " * id()) + 2] *= 2.0;\n"
		<< "    totalVals[(" << 2 * nVars << " * id()) + 5] *= 2.0;\n"
		<< "  }\n"
		   "}\n"
		   "\n";

	ofs << "std::size_t\n"
		<< _reactionNetwork
		<< "::getMonitorDataLineSize() const\n"
		   "{\n"
		<< "  return getSpeciesListSize() * " << 2 * nVars << ";\n"
		<< "}\n"
		   "\n";

	ofs << "void\n"
		<< _reactionNetwork << "::writeMonitorDataLine(\n"
		<< "  const std::vector<double>& localData, double time)\n"
		   "{\n"
		   "  auto numSpecies = getSpeciesListSize();\n"
		   "  auto globalData = std::vector<double>(localData.size(), 0.0);\n"
		   "  MPI_Reduce(localData.data(), globalData.data(),\n"
		   "    localData.size(), MPI_DOUBLE, MPI_SUM, 0,\n"
		   "    util::getMPIComm());\n"
		   "\n"
		   "  if (util::getMPIRank() == 0) {\n"
		   "    for (auto i = 0; i < numSpecies; ++i) {\n"
		<< "      auto id = [i](std::size_t n)\n"
		<< "        { return " << 2 * nVars << " * i + n; };\n"
		<< "      constexpr double tol = 1.0e-16;\n"
		<< "      if (globalData[id(0)] > tol) {\n"
		<< "        globalData[id(" << nVars - 1
		<< ")] /= globalData[id(0)];\n"
		   "      }\n"
		<< "      if (globalData[id(" << nVars << ")] > tol) {\n"
		<< "        globalData[id(" << 2 * nVars - 1 << ")] /=\n"
		<< "          globalData[id(" << nVars << ")];\n"
		<< "      }\n"
		   "    }\n"
		   "    const int outputPrecision = 5;\n"
		   "    std::fstream outputFile;\n"
		   "    outputFile.open(getMonitorOutputFileName(),\n"
		   "      std::fstream::out | std::fstream::app);\n"
		   "    outputFile << std::setprecision(outputPrecision);\n"
		   "    outputFile << time;\n"
		   "    for (auto i = 0; i < numSpecies; ++i) {\n"
		   "      auto id = [i](std::size_t n)\n"
		<< "        { return " << 2 * nVars << " * i + n; };\n";
	for (auto di = 0; di < nVars; ++di) {
		ofs << "      outputFile << ' ' << globalData[id(" << di << ")];\n";
	}
	for (auto di = 0; di < nVars; ++di) {
		ofs << "      outputFile << ' ' << globalData[id(" << nVars + di
			<< ")];\n";
	}
	ofs << "    }\n"
		   "    outputFile << std::endl;\n"
		   "    outputFile.close();\n"
		   "  }\n"
		   "}\n";

	ofs << "}\n";
}

void
NetworkHandlerClassGenerator::generateNetwork()
{
	generateNetworkHeader();
	generateNetworkImpl();
}

void
NetworkHandlerClassGenerator::generateNetworkHandler()
{
	auto filePath = _genDir / (_networkHandler + ".h");
	auto ofs = openFile(filePath);
	ofs << "#pragma once\n"
		   "#include <xolotl/core/network/NetworkHandler.h>\n"
		   "#include <xolotl/factory/network/NetworkHandlerFactory.h>\n"
		   "\n"
		   "namespace xolotl::core::network\n"
		   "{\n"
		<< "class " << _networkHandler << " : public NetworkHandler\n"
		<< "{\n"
		<< "public:\n"
		<< "  " << _networkHandler << "(const options::IOptions& options);\n"
		<< "};\n"
		<< "}\n";

	filePath = _genDir / (_networkHandler + ".cpp");
	ofs = openFile(filePath);
	ofs << "#include <xolotl/factory/network/NetworkHandlerFactory.h>\n"
		<< "#include <" << _networkHandler << ".h>\n"
		<< "#include <" << _reactionNetwork << ".h>\n"
		<< "\n"
		   "namespace xolotl::core::network\n"
		   "{\n"
		   "namespace detail\n"
		   "{\n"
		   "using NHF = ::xolotl::factory::network::NetworkHandlerFactory;\n"
		   "template <typename T>\n"
		   "using RegCol = NHF::RegistrationCollection<T>;\n"
		   "auto networkHandlerRegistrations =\n"
		<< "  RegCol<" << _networkHandler << ">(\n"
		<< "    {\"" << _userMaterial << "\"});\n"
		<< "}\n"
		   "\n"
		   "auto nwGenerator = [](const options::IOptions& options) {\n"
		<< "  using NetworkType = " << _reactionNetwork << ";\n"
		<< "  using AmountType = NetworkType::AmountType;\n"
		   "  AmountType maxV = options.getMaxV();\n"
		   "  AmountType maxI = options.getMaxI();\n"
		   "  std::vector<AmountType> maxSpeciesAmounts = {\n"
		   "    maxV, maxI};\n"
		   "  std::vector<NetworkType::SubdivisionRatio> subdivRatios = {\n"
		   "    {maxV + 1, maxI + 1}};\n"
		   "  auto network = std::make_shared<NetworkType>(\n"
		   "    maxSpeciesAmounts, subdivRatios, 1, options);\n"
		   "  return network;\n"
		   "};\n"
		   "\n"
		<< _networkHandler << "::" << _networkHandler << "(\n"
		<< "  const options::IOptions& options)\n"
		<< "  : NetworkHandler(options, nwGenerator)\n"
		<< "{\n"
		<< "}\n";

	ofs << "}\n";

	ofs << "#include <xolotl/factory/material/MaterialHandlerFactory.h>\n"
		<< "#include <xolotl/core/flux/CustomFitFluxHandler.h>\n"
		<< "#include <xolotl/core/material/MaterialHandler.h>\n";

	ofs << "namespace xolotl::core::material\n"
		<< "{\n"
		<< "class " << _materialHandler << " : public MaterialHandler\n"
		<< "{\n"
		<< "public:\n"
		<< "  using SubHandlerGenerator =\n"
		<< "    MaterialSubHandlerGenerator<flux::CustomFitFluxHandler>;\n"
		<< "  " << _materialHandler << "(const options::IOptions& options) :\n"
		<< "    MaterialHandler(options, SubHandlerGenerator{})\n"
		<< "  {}\n"
		<< "};\n"
		<< "\n"
		<< "namespace detail\n"
		   "{\n"
		   "using MHF = ::xolotl::factory::material::MaterialHandlerFactory;\n"
		   "template <typename T>\n"
		   "using RegCol = MHF::RegistrationCollection<T>;\n"
		   "auto materialHandlerRegistrations = RegCol<"
		<< _materialHandler << ">(\n"
		<< "    {\"" << _userMaterial << "\"});\n"
		<< "}\n"
		<< "}\n";
}

void
NetworkHandlerClassGenerator::writeCMakeLists()
{
	auto filePath = _genDir / "CMakeLists.txt";
	auto ofs = openFile(filePath);

	ofs << "cmake_minimum_required(VERSION 3.13)\n"
		<< "project(" << _networkHandler << " LANGUAGES CXX)\n"
		<< "set(CMAKE_CXX_STANDARD 17)\n"
#ifndef NDEBUG
		<< "set(CMAKE_BUILD_TYPE Debug)\n"
#endif
		<< "set(CMAKE_VERBOSE_MAKEFILE ON)\n"
		<< "find_package(Xolotl)\n"
		<< "add_library(" << _networkHandler << " SHARED\n"
		<< _networkHandler << ".cpp\n"
		<< _reactionNetwork << ".cpp\n"
		<< _clusterGenerator << ".cpp\n"
		<< ")\n"
		<< "target_link_libraries(" << _networkHandler << " xolotlCore)\n"
		<< "target_include_directories(" << _networkHandler
		<< " PUBLIC ${CMAKE_CURRENT_LIST_DIR})\n";
}

void
NetworkHandlerClassGenerator::generateBuild()
{
	writeCMakeLists();

	auto execDir = _execFile.parent_path();
	auto xolotlPrefix = execDir.parent_path();
	fs::create_directories(_buildDir);
	XOLOTL_LOG_XTRA << "XNGEN: xolotl install prefix: " << xolotlPrefix;
	std::stringstream cmdss;
	cmdss << "cmake -S " << _genDir << " -B " << _buildDir
		  << " -DCMAKE_PREFIX_PATH=" << xolotlPrefix;

	auto outputFile = _buildDir / "out.txt";
	fs::remove(outputFile);
	auto redirect = " >> " + outputFile.string() + " 2>&1";

	XOLOTL_LOG_XTRA << "XNGEN: config command: " << cmdss.str();
	std::cout << std::flush;
	auto rv = std::system((cmdss.str() + redirect).c_str());
	if (rv != 0) {
		XOLOTL_ERROR(std::runtime_error,
			"config step failed; see " + outputFile.string());
	}
	XOLOTL_LOG_XTRA << "XNGEN: config return: " << rv;
	cmdss.str("");
	auto nthreads = std::thread::hardware_concurrency();
	nthreads = std::lround(std::ceil(0.25 * nthreads));
	cmdss << "cmake --build " << _buildDir << " -j " << nthreads;
	XOLOTL_LOG_XTRA << "XNGEN: build command: " << cmdss.str();
	std::cout << std::flush;
	rv = std::system((cmdss.str() + redirect).c_str());
	if (rv != 0) {
		XOLOTL_ERROR(std::runtime_error,
			"build step failed; see " + outputFile.string());
	}
	XOLOTL_LOG_XTRA << "XNGEN: build return: " << rv;
}

void
loadLibrary(const std::string& libFile)
{
	boost::dll::shared_library lib(libFile);
}

void
broadcastString(std::string& str)
{
	auto rank = util::getMPIRank();
	auto comm = util::getMPIComm();
	int len = (rank == 0) ? static_cast<int>(str.length()) : -1;
	MPI_Bcast(&len, 1, MPI_INT, 0, comm);
	if (rank != 0) {
		str = std::string(len, ' ');
	}
	MPI_Bcast(str.data(), len + 1, MPI_CHAR, 0, comm);
}

bool
NetworkHandlerClassGenerator::needToGenerate() const
{
	if (!exists(_networkLibFile)) {
		return true;
	}
	auto rnFileMTime = last_write_time(_rnFile);
	auto exeBuildTime = last_write_time(_execFile);
	auto libBuildTime = last_write_time(_networkLibFile);
	return (exeBuildTime > libBuildTime) || (libBuildTime < rnFileMTime);
}

template <typename F>
void
time(F&& f, const std::string& message)
{
	XOLOTL_LOG << message << "... ";
	auto start = std::chrono::steady_clock::now();
	f();
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	XOLOTL_LOG << "  " << elapsed_seconds.count() << "s\n";
	::xolotl::util::Log::flush();
}

void
NetworkHandlerClassGenerator::generate()
{
	if (!this->needToGenerate()) {
		return;
	}
	time(
		[this] {
			generateTraits();
			generateReactions();
			generateClusterGenerator();
			generateNetwork();
			generateNetworkHandler();
		},
		"Generate custom network code");

	time([this] { this->generateBuild(); }, "Build custom network");
}

void
generateNetworkHandler(options::IOptions* options)
{
	std::string libFile;

	if (util::getMPIRank() == 0) {
		const auto& rnFile = options->getReactionNetworkFileName();
		const auto& material = options->getMaterial();
		auto userMat = material + "_u";
		options->setMaterial(userMat);
		auto gen = NetworkHandlerClassGenerator{rnFile, material, userMat};
		libFile = gen.getLibraryFileName();
		gen.generate();
	}

	broadcastString(libFile);

	if (util::getMPIRank() == 0) {
		time([&] { loadLibrary(libFile); }, "Load custom network");
	}
	else {
		loadLibrary(libFile);
	}
}
} // namespace interface
} // namespace xolotl
