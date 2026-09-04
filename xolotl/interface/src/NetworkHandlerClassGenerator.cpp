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
#include <boost/rational.hpp>

#include <xolotl/interface/NetworkHandlerClassGenerator.h>
#include <xolotl/util/Filesystem.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/StreamUtils.h>
#include <xolotl/util/Tokenizer.h>

using namespace std::string_literals;

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

auto xOr = [](auto&& a, auto&& b) { return !a != !b; };

NetworkHandlerClassGenerator::NetworkHandlerClassGenerator(
	const fs::path& rnFile, const std::string& material,
	const std::string& userMaterial) :
	_rnFile(rnFile),
	_material(material),
	_userMaterial(userMaterial),
	_genDir(fs::current_path() / ".xgrn"),
	_buildDir(_genDir / "build"),
	_execFile(boost::dll::program_location().string())
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

template <typename TArray, typename TEnum>
struct EnumIndexed : TArray
{
	using DimType = std::underlying_type_t<TEnum>;
	decltype(auto)
	operator[](TEnum enumVal)
	{
		return TArray::operator[](static_cast<DimType>(enumVal));
	}
};

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

	/*
	 * network
	 */
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

	using SpType = SpeciesData::Type;
	EnumIndexed<std::array<std::vector<SpeciesData>, 3>, SpType> spTypeMap;
	auto mapSpeciesType = [&](const std::string& type) -> SpeciesData& {
		if (type == "impurity") {
			auto& sp = spTypeMap[SpType::impurity].emplace_back();
			sp.type = SpType::impurity;
			return sp;
		}
		else if (type == "vacancy") {
			auto& sp = spTypeMap[SpType::vacancy].emplace_back();
			sp.type = SpType::vacancy;
			return sp;
		}
		else if (type == "interstitial") {
			auto& sp = spTypeMap[SpType::interstitial].emplace_back();
			sp.type = SpType::interstitial;
			return sp;
		}
		else {
			XOLOTL_ERROR(std::invalid_argument,
				"unsupported species type: \"" + type + "\"");
		}
	};

	/*
	 * species
	 */
	auto spNode = tree.get_child("species");
	for (auto&& [_, elNode] : spNode) {
		auto& spData = mapSpeciesType(toLower(elNode.get<std::string>("type")));
		spData.name = elNode.get<std::string>("name");
		spData.label = elNode.get<std::string>("label");
		if (elNode.count("max_pure_amount")) {
			spData.maxPure = elNode.get<AmountType>("max_pure_amount");
		}
		if (elNode.count("interact")) {
			auto ixNode = elNode.get_child("interact");
			auto ixSpecies = ixNode.get<std::string>("species");
			auto ixR = ixNode.get<boost::rational<AmountType>>("ratio");
			spData.interact.push_back({.species = ixSpecies, .ratio = ixR});
		}
	}
	const auto& impSp = spTypeMap[SpType::impurity];
	std::copy(begin(impSp), end(impSp), std::back_inserter(_speciesData));
	const auto& vacSp = spTypeMap[SpType::vacancy];
	std::copy(begin(vacSp), end(vacSp), std::back_inserter(_speciesData));
	const auto& intSp = spTypeMap[SpType::interstitial];
	std::copy(begin(intSp), end(intSp), std::back_inserter(_speciesData));
	for (auto&& sp : _speciesData) {
		if (!_speciesLabelMap.try_emplace(sp.label, &sp).second) {
			XOLOTL_ERROR(std::invalid_argument,
				"Duplicate species label \"" + sp.label + "\"");
		}
	}

	auto checkSpeciesLabel = [*this](const std::string& label,
								 const std::string& msgPrefix) {
		if (!_speciesLabelMap.count(label)) {
			XOLOTL_ERROR(std::invalid_argument,
				msgPrefix + " species (" + label +
					") not found in species list");
		}
	};

	/*
	 ************
	 * clusters *
	 ************
	 */

	// parse label string as one species or two separated with a hyphen
	auto parseClusterType = [&](const std::string& label) {
		auto pos = label.find("-");
		if (pos == std::string::npos) {
			checkSpeciesLabel(label, "cluster");
			return std::array{label, ""s};
		}
		auto label0 = label.substr(0, pos);
		auto label1 = label.substr(pos + 1);
		checkSpeciesLabel(label0, "cluster");
		checkSpeciesLabel(label1, "cluster");
		return std::array{label0, label1};
	};

	auto getClusterConstituents = [this](
									  const std::array<std::string, 2> labels) {
		const SpeciesData* sp0 = _speciesLabelMap.at(labels[0]);
		const SpeciesData* sp1 = nullptr;
		if (!labels[1].empty()) {
			sp1 = _speciesLabelMap.at(labels[1]);
		}
		return std::array{sp0, sp1};
	};

	// read clusters
	auto clNode = tree.get_child("clusters");
	for (auto&& [_, elNode] : clNode) {
		auto& clData = _clusterData.emplace_back();
		// cluster type
		auto typeLabel = elNode.get<std::string>("type");
		clData.type.label = typeLabel;
		clData.type.constituentLabels = parseClusterType(typeLabel);
		clData.type.constituents =
			getClusterConstituents(clData.type.constituentLabels);
		if (clData.type.constituents[1] == nullptr) {
			clData.type.species = clData.type.constituents[0];
		}
		// cluster size range
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
		// cluster radius
		if (!elNode.count("radius")) {
			XOLOTL_ERROR(std::invalid_argument,
				"every cluster record must specify \"radius\"");
		}
		clData.radiusExpr = elNode.get<std::string>("radius");
		// cluster parameters
		clData.migrationEnergy = elNode.get("migration_energy", nan);
		clData.diffusionFactor = elNode.get("diffusion_factor", nan);
	}
	// group clusters by type
	for (const auto& cl : _clusterData) {
		auto& grp = _clusterGroups[cl.type.label];
		if (grp.type.label.empty()) {
			grp.type = cl.type;
		}
		if (cl.size[0] == maxAmount && cl.size[1] == maxAmount) {
			if (grp.general != nullptr) {
				XOLOTL_ERROR(std::runtime_error,
					"multiple definitions of general case for " +
						cl.type.label + " type clusters");
			}
			grp.general = &cl;
		}
		else {
			grp.cases.push_back(&cl);
		}
	}

	auto checkClusterTypeLabel = [*this](const std::string& label,
									 const std::string& msgPrefix) {
		if (!_clusterGroups.count(label) && label != "0") {
			XOLOTL_ERROR(std::invalid_argument,
				msgPrefix + " cluster type (" + label + ") not found");
		}
	};

	/*
	 *************
	 * reactions *
	 *************
	 */

	// parse expression to get reactant and product labels
	auto parseReactionExpr = [checkClusterTypeLabel](const std::string& expr) {
		auto lr = util::Tokenizer<>{expr, "->"}();
		auto reactants = util::Tokenizer<>{lr[0], "+"}();
		auto products = util::Tokenizer<>{lr[1], "+"}();
		for (auto& r : reactants) {
			r = removeSpaces(r);
			checkClusterTypeLabel(r, "reactant");
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
			checkClusterTypeLabel(p, "product");
		}
		return std::make_tuple(reactants, products);
	};

	// read reactions
	auto rnNode = tree.get_child("reactions");
	for (auto&& [_, elNode] : rnNode) {
		auto& rnData = _reactionData.emplace_back();
		if (elNode.count("expr")) {
			std::tie(rnData.parts.reactantLabels, rnData.parts.productLabels) =
				parseReactionExpr(elNode.get<std::string>("expr"));
		}
		else {
			if (!(elNode.count("reactants") && elNode.count("products"))) {
				XOLOTL_ERROR(std::runtime_error,
					"Reaction must be expressed either as \"expr\" or as "
					"\"reactants\" and \"products\"");
			}

			// get reactant labels from list
			auto reactantNode = elNode.get_child("reactants");
			if (reactantNode.empty()) {
				XOLOTL_ERROR(std::runtime_error, "Reactant list empty");
			}
			for (auto&& reac : reactantNode) {
				auto s = reac.second.get_value<std::string>();
				checkClusterTypeLabel(s, "reactant");
				rnData.parts.reactantLabels.push_back(s);
			}

			// get product labels from list
			auto productNode = elNode.get_child("products");
			if (productNode.empty()) {
				XOLOTL_ERROR(std::runtime_error, "Product list empty");
			}
			for (auto&& prod : productNode) {
				auto s = prod.second.get_value<std::string>();
				checkClusterTypeLabel(s, "product");
				rnData.parts.productLabels.push_back(s);
			}
		}
		for (auto&& r : rnData.parts.reactantLabels) {
			rnData.parts.reactants.push_back(_clusterGroups.at(r).type);
		}
		if (rnData.parts.productLabels == std::vector{"0"s}) { }
		else {
			for (auto&& p : rnData.parts.productLabels) {
				rnData.parts.products.push_back(_clusterGroups.at(p).type);
			}
		}
		rnData.size = elNode.get("size", maxAmount);
		if (elNode.count("binding")) {
			rnData.bindingExpr = elNode.get<std::string>("binding");
		}
	}
	// group reactions by species set
	for (const auto& rn : _reactionData) {
		if (rn.parts.reactants.size() > 1) {
			std::string key{};
			for (const auto& r : rn.parts.reactantLabels) {
				key += r;
			}
			for (const auto& p : rn.parts.productLabels) {
				key += p;
			}
			auto& grp = _productionReactionGroups[key];
			grp.parts = rn.parts;
			if (rn.size == maxAmount) {
				grp.general = &rn;
			}
			else {
				grp.cases.push_back(&rn);
			}
		}
		else {
			std::string key = rn.parts.reactantLabels[0];
			for (const auto& p : rn.parts.productLabels) {
				key += p;
			}
			auto& grp = _dissociationReactionGroups[key];
			grp.parts = rn.parts;
			if (rn.size == maxAmount) {
				grp.general = &rn;
			}
			else {
				grp.cases.push_back(&rn);
			}
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
	for (auto&& spec : _speciesData) {
		for (auto&& ix : spec.interact) {
			auto spLab = spec.label;
			auto ispLab = ix.species;
			auto ixLab = spLab + ispLab;
			ofs << "  AmountType _ix_" << ixLab << "_max_" << spLab << ";\n"
				<< "  AmountType _ix_" << ixLab << "_max_" << ispLab << ";\n";
		}
	}
	// ofs << "  AmountType _maxV{0};\n";
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
		   "#include <ratio>\n"
		   "\n"
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
		   "  for (auto& res : result) {\n"
		   "    res = true;\n"
		   "  }\n"
		   "  return true;\n"
		   "}\n"
		   "\n"
		   "KOKKOS_INLINE_FUNCTION\n"
		   "bool\n"
		<< _clusterGenerator << "::select(const Region& region) const\n"
		<< "{\n"
		   "  constexpr auto species = NetworkType::getSpeciesRange();\n"
		   "  auto matchSpeciesPair =\n"
		   "    [&species](Composition lo, Species s1, Species s2) {\n"
		   "      if (lo[s1] == 0 || lo[s2] == 0) { return false; }\n"
		   "      for (auto l : species) {\n"
		   "        if (l == s1 || l == s2) { continue; }\n"
		   "        if (lo[l] != 0) { return false; }\n"
		   "      }\n"
		   "      return true;\n"
		   "    };\n"
		   "\n"
		   "  auto origin = Composition(region.getOrigin());\n";

	for (auto&& spec : _speciesData) {
		auto spRef = "Species::" + spec.label;
		auto mxp = spec.maxPure;
		ofs << "  if (region.isSimplex() &&\n"
			<< "      origin.isOnAxis(" << spRef << ")) {\n";

		if (mxp != bad<AmountType>) {
			ofs << "    if (origin[" << spRef << "] <= " << mxp << ") {\n"
				<< "      return true;\n"
				<< "    }\n";
		}
		else {
			ofs << "    return true;\n";
		}
		ofs << "  }\n";
		for (auto&& ix : spec.interact) {
			auto ispRef = "Species::" + ix.species;
			auto ixLab = spec.label + ix.species;
			auto ixSpMax = "_ix_" + ixLab + "_max_" + spec.label;
			auto ixISpMax = "_ix_" + ixLab + "_max_" + ix.species;
			ofs << "  {\n"
				<< "  constexpr auto sp = " << spRef << ";\n"
				<< "  constexpr auto ix = " << ispRef << ";\n"
				<< "  if (matchSpeciesPair(origin, sp, ix)) {\n"
				<< "    if (origin[sp] <= " << ixSpMax << " &&\n"
				<< "        origin[ix] <= " << ixISpMax << ") {\n"
				<< "      return true;\n"
				<< "    }\n"
				<< "    return false;\n"
				<< "  }\n"
				<< "  }\n";
		}
	}
	ofs << "  return false;\n"
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
	for (auto&& [type, grp] : _clusterGroups) {
		if (grp.cases.empty()) {
			continue;
		}
		if (grp.type.species == nullptr) {
			// TODO
			continue;
		}
		auto species = "Species::" + type;
		ofs << "    if (comp.isOnAxis(" << species << ")) {\n"
			<< "      auto amt = comp[" << species << "];\n";
		for (auto cl : grp.cases) {
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
		}
		if (grp.general) {
			if (!std::isnan(grp.general->migrationEnergy)) {
				ofs << "      return " << grp.general->migrationEnergy << ";\n";
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
	for (auto&& [type, grp] : _clusterGroups) {
		if (grp.cases.empty()) {
			continue;
		}
		if (grp.type.species == nullptr) {
			// TODO
			continue;
		}
		auto species = "Species::" + type;
		ofs << "    if (comp.isOnAxis(" << species << ")) {\n"
			<< "      auto amt = comp[" << species << "];\n";
		for (auto cl : grp.cases) {
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
		}
		if (grp.general) {
			if (!std::isnan(grp.general->diffusionFactor)) {
				ofs << "      return " << grp.general->diffusionFactor << ";\n";
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
	for (auto&& [_, grp] : _clusterGroups) {
		if (grp.cases.empty()) {
			// TODO: do I need this?
			continue;
		}
		if (grp.type.species == nullptr) {
			auto s1 = grp.type.constituentLabels[0];
			auto s2 = grp.type.constituentLabels[1];
			auto spec1 = "Species::" + s1;
			auto spec2 = "Species::" + s2;
			auto amt1 = "amt" + s1;
			auto amt2 = "amt" + s2;
			ofs << "    {\n"
				<< "    auto " << amt1 << " = comp[" << spec1 << "];\n"
				<< "    auto " << amt2 << " = comp[" << spec2 << "];\n"
				<< "    if (\n"
				<< "        " << amt1 << " > 0\n"
				<< "        && " << amt2 << " > 0";
			for (const auto& sp : _speciesData) {
				if (sp.label == s1 || sp.label == s2) {
					continue;
				}
				ofs << "        && comp[Species::" << sp.label << "] == 0\n";
			}
			ofs << "        ) {\n";
			for (auto cl : grp.cases) {
				const auto& rExpr = cl->radiusExpr;
				auto expr = parseRadiusExpr(rExpr);
				ofs << "      return " << expr << ";\n";
			}
			ofs << "    }\n"
				<< "    }\n";
		}
		else {
			auto species = "Species::" + grp.type.label;
			ofs << "    if (comp.isOnAxis(" << species << ")) {\n"
				<< "      auto amt = comp[" << species << "];\n";
			for (auto cl : grp.cases) {
				const auto& rExpr = cl->radiusExpr;
				auto expr = parseRadiusExpr(rExpr);
				if (cl->size[0] < maxAmount) {
					if (cl->size[1] < maxAmount) {
						ofs << "      if (" << cl->size[0] << " <= amt &&\n"
							<< "          amt <= " << cl->size[1] << ") {\n"
							<< "        // parsed from: " << rExpr << "\n"
							<< "        return " << expr << ";\n"
							<< "      }\n";
					}
					else {
						ofs << "      if (amt == " << cl->size[0] << ") {\n"
							<< "        // parsed from: " << rExpr << "\n"
							<< "        return " << expr << ";\n"
							<< "      }\n";
					}
				}
			}
			if (grp.general) {
				const auto& rExpr = grp.general->radiusExpr;
				auto expr = parseRadiusExpr(rExpr);
				ofs << "      // parsed from: " << rExpr << "\n"
					<< "      return " << expr << ";\n";
			}
			ofs << "    }\n";
		}
	}
	ofs << "  }\n"
		   "  return 0.0;\n"
		   "}\n"
		   "}\n";

	filePath = _genDir / (_clusterGenerator + ".cpp");
	ofs = openFile(filePath);
	ofs << "#include <ratio>\n"
		<< "#include <" << _reactionNetwork << ".h>\n"
		<< "namespace xolotl::core::network {\n"
		<< _clusterGenerator << "::" << _clusterGenerator << "(\n"
		<< "  const options::IOptions& options)\n"
		<< "{\n"
		   "  const auto& params = options.getNetworkParameters();\n";
	for (auto&& spec : _speciesData) {
		for (auto&& ix : spec.interact) {
			auto spLab = spec.label;
			auto spRef = "Species::" + spLab;
			auto ispLab = ix.species;
			auto ispRef = "Species::" + ispLab;
			auto ixLab = spLab + ispLab;
			auto spMaxVar = "_ix_" + ixLab + "_max_" + spLab;
			auto spMaxOpt = "max" + spLab;
			auto ispMaxVar = "_ix_" + ixLab + "_max_" + ispLab;
			auto ispMaxOpt = "max" + ispLab;
			ofs << "  {\n"
				<< "  AmountType " << spMaxOpt << " =\n";
			switch (spec.type) {
			case SpeciesData::Type::impurity:
				ofs << "    params[0];\n";
				break;
			case SpeciesData::Type::vacancy:
				ofs << "    params[3];\n";
				break;
			case SpeciesData::Type::interstitial:
				ofs << "    params[4];\n";
				break;
			}
			ofs << "  AmountType " << ispMaxOpt << " =\n";
			switch (_speciesLabelMap[ispLab]->type) {
			case SpeciesData::Type::impurity:
				ofs << "    params[0];\n";
				break;
			case SpeciesData::Type::vacancy:
				ofs << "    params[3];\n";
				break;
			case SpeciesData::Type::interstitial:
				ofs << "    params[4];\n";
				break;
			}
			ofs << "  constexpr auto iR = std::ratio<\n"
				<< "    " << ix.ratio.numerator() << ",\n"
				<< "    " << ix.ratio.denominator() << ">{};\n"
				<< "  if (iR.num <= iR.den) {\n"
				<< "    " << ispMaxVar << " = " << ispMaxOpt << ";\n"
				<< "    " << spMaxVar << " = " << ispMaxOpt << " *\n"
				<< "      iR.num / iR.den;\n"
				<< "    " << spMaxVar << " = std::max(\n"
				<< "      " << spMaxVar << ", " << spMaxOpt << ");\n"
				<< "  }\n"
				<< "  else {\n"
				<< "    " << spMaxVar << " = " << spMaxOpt << ";\n"
				<< "    " << ispMaxVar << " = " << spMaxOpt << " *\n"
				<< "      iR.den / iR.num;\n"
				<< "    " << ispMaxVar << " = std::max(\n"
				<< "      " << ispMaxVar << ", " << ispMaxOpt << ");\n"
				<< "  }\n"
				<< "  }\n";
		}
	}
	ofs << "}\n"
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
		   "  getRate(const TRegion& pairCl0Reg, const TRegion& "
		   "pairCl1Reg, \n"
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

	auto expressCases = [&](auto&& rnGroup) {
		std::stringstream oss;
		for (auto&& rnCase : rnGroup.cases) {
			const auto& expr = rnCase->bindingExpr;
			oss << "    if (amt == " << rnCase->size << ") {\n"
				<< "      be = " << parseBindingExpr(expr) << ";\n"
				<< "    }\n"
				<< "    else ";
		}
		oss << "    {\n";
		if (rnGroup.general) {
			const auto& expr = rnGroup.general->bindingExpr;
			oss << "      be = " << parseBindingExpr(expr) << ";\n";
		}
		oss << "    }\n";
		return oss.str();
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
		   "      prod2Reg.isSimplex()) {\n"
		   "    Composition rComp = clReg.getOrigin();\n"
		   "    Composition p1Comp = prod1Reg.getOrigin();\n"
		   "    Composition p2Comp = prod2Reg.getOrigin();\n"
		   "\n"
		   "  constexpr auto species = NetworkType::getSpeciesRange();\n"
		   "  auto matchSpeciesPair =\n"
		   "    [&species](Composition comp, Species s1, Species s2) {\n"
		   "      if (comp[s1] == 0 || comp[s2] == 0) { return false; }\n"
		   "      for (auto l : species) {\n"
		   "        if (l == s1 || l == s2) { continue; }\n"
		   "        if (comp[l] != 0) { return false; }\n"
		   "      }\n"
		   "      return true;\n"
		   "    };\n";
	for (auto&& [key, rnGroup] : _dissociationReactionGroups) {
		const auto& reactants = rnGroup.parts.reactants;
		auto r = reactants[0];
		auto p1 = rnGroup.parts.products[0];
		auto p2 = rnGroup.parts.products[1];
		if (reactants.size() == 1 && r.species == nullptr) {
			if (!rnGroup.general) {
				XOLOTL_ERROR(std::invalid_argument,
					"no general case provided for reaction " + key);
			}
			const auto& expr = rnGroup.general->bindingExpr;
			auto rSp1Label = r.constituentLabels[0];
			auto rSp2Label = r.constituentLabels[1];
			auto rSp1Ref = "Species::" + rSp1Label;
			auto rSp2Ref = "Species::" + rSp2Label;

			if (p1.species && p2.species) {
				auto p1SpRef = "Species::" + p1.label;
				auto p2SpRef = "Species::" + p2.label;
				ofs << "    if (matchSpeciesPair(rComp,\n"
					<< "          " << rSp1Ref << ", " << rSp2Ref << ") &&\n"
					<< "        (p1Comp.isOnAxis(" << p1SpRef << ") &&\n"
					<< "        p2Comp.isOnAxis(" << p2SpRef << ")) ||\n"
					<< "        (p1Comp.isOnAxis(" << p2SpRef << ") &&\n"
					<< "        p2Comp.isOnAxis(" << p1SpRef << "))) {\n"
					<< "      auto amt" << rSp1Label << " =\n"
					<< "        rComp[" << rSp1Ref << "];\n"
					<< "      auto amt" << rSp2Label << " =\n"
					<< "        rComp[" << rSp2Ref << "];\n"
					<< "      // parsed from: " << expr << "\n"
					<< "      be = " << parseBindingExpr(expr) << ";\n"
					<< "    }\n";
			}
			else if (xOr(p1.species, p2.species)) {
				std::string pfx = "Species::";
				auto [p1SpRef, p2Sp1Ref, p2Sp2Ref] = p1.species ?
					std::make_tuple(pfx + p1.label,
						pfx + p2.constituentLabels[0],
						pfx + p2.constituentLabels[1]) :
					std::make_tuple(pfx + p2.label,
						pfx + p1.constituentLabels[0],
						pfx + p1.constituentLabels[1]);
				ofs << "    if (matchSpeciesPair(rComp,\n"
					<< "          " << rSp1Ref << ", " << rSp2Ref << ") &&\n"
					<< "        (p1Comp.isOnAxis(" << p1SpRef << ") &&\n"
					<< "        matchSpeciesPair(p2Comp,\n"
					<< "          " << p2Sp1Ref << ", " << p2Sp2Ref << ")) ||\n"
					<< "        (p2Comp.isOnAxis(" << p1SpRef << ") &&\n"
					<< "        matchSpeciesPair(p1Comp,\n"
					<< "          " << p2Sp1Ref << ", " << p2Sp2Ref << "))) {\n"
					<< "      auto amt" << rSp1Label << " =\n"
					<< "        rComp[" << rSp1Ref << "];\n"
					<< "      auto amt" << rSp2Label << " =\n"
					<< "        rComp[" << rSp2Ref << "];\n"
					<< "      // parsed from: " << expr << "\n"
					<< "      be = " << parseBindingExpr(expr) << ";\n"
					<< "    }\n";
			}
			else {
				auto p1Sp1Ref = "Species::" + p1.constituentLabels[0];
				auto p1Sp2Ref = "Species::" + p1.constituentLabels[1];
				auto p2Sp1Ref = "Species::" + p2.constituentLabels[0];
				auto p2Sp2Ref = "Species::" + p2.constituentLabels[1];
				ofs << "    if (matchSpeciesPair(rComp,\n"
					<< "          " << rSp1Ref << ", " << rSp2Ref << ") &&\n"
					<< "        (matchSpeciesPair(p1Comp,\n"
					<< "          " << p2Sp1Ref << ", " << p2Sp2Ref << ") &&\n"
					<< "        matchSpeciesPair(p2Comp,\n"
					<< "          " << p1Sp1Ref << ", " << p1Sp2Ref << ")) ||\n"
					<< "        (matchSpeciesPair(p2Comp,\n"
					<< "          " << p2Sp1Ref << ", " << p2Sp2Ref << ") &&\n"
					<< "        matchSpeciesPair(p1Comp,\n"
					<< "          " << p1Sp1Ref << ", " << p1Sp2Ref << "))) {\n"
					<< "      auto amt" << rSp1Label << " =\n"
					<< "        rComp[" << rSp1Ref << "];\n"
					<< "      auto amt" << rSp2Label << " =\n"
					<< "        rComp[" << rSp2Ref << "];\n"
					<< "      // parsed from: " << expr << "\n"
					<< "      be = " << parseBindingExpr(expr) << ";\n"
					<< "    }\n";
			}
			ofs << "    else";
			continue;
		}
		if (reactants.size() == 1 && xOr(p1.species, p2.species)) {
			auto rSpLabel = r.label;
			std::string pfx = "Species::";
			auto rSpRef = pfx + reactants[0].label;
			auto [p1SpRef, p2Sp1Ref, p2Sp2Ref] = p1.species ?
				std::make_tuple(pfx + p1.label, pfx + p2.constituentLabels[0],
					pfx + p2.constituentLabels[1]) :
				std::make_tuple(pfx + p2.label, pfx + p1.constituentLabels[0],
					pfx + p1.constituentLabels[1]);
			ofs << "    if (rComp.isOnAxis(" << rSpRef << ") &&\n"
				<< "        (p1Comp.isOnAxis(" << p1SpRef << ") &&\n"
				<< "          matchSpeciesPair(p2Comp,\n"
				<< "            " << p2Sp1Ref << ", " << p2Sp2Ref << ")) ||\n"
				<< "        (p2Comp.isOnAxis(" << p1SpRef << ") &&\n"
				<< "          matchSpeciesPair(p1Comp,\n"
				<< "            " << p2Sp2Ref << ", " << p2Sp2Ref << "))) {\n"
				<< "      auto amt = rComp[" << rSpRef << "];\n"
				<< expressCases(rnGroup)
				// << "      be = " << parseBindingExpr(expr) << ";\n"
				<< "    }\n";
			ofs << "    else";
			continue;
		}
		if (reactants.size() == 1 && p1.label == p2.label) {
			auto spec = p1.label;
			auto specRef = "Species::" + spec;
			auto rSpLabel = r.label;
			auto rSpRef = "Species::" + rSpLabel;
			ofs << "    if (rComp.isOnAxis(" << rSpRef << ") &&\n"
				<< "        p1Comp.isOnAxis(" << specRef << ") &&\n"
				<< "        p2Comp.isOnAxis(" << specRef << ")) {\n"
				<< "      auto amt = rComp[" << rSpRef << "];\n"
				<< expressCases(rnGroup)
				// << "      be = " << parseBindingExpr(expr) << ";\n"
				<< "    }\n";
			ofs << "    else";
			continue;
		}
		XOLOTL_ERROR(std::runtime_error, "unsupported reaction case");
	}
	ofs << "    {}\n"
		   "  }\n"
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
		   "  KOKKOS_INLINE_FUNCTION\n"
		   "  void\n"
		   "  setConnectivity(Connectivity)\n"
		   "  {\n"
		   "  }\n"
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
		   "  constexpr auto speciesNoI = "
		   "    NetworkType::getSpeciesRangeNoI();\n"
		   "\n"
		   "  const auto& cl1Reg = this->getCluster(i).getRegion();\n"
		   "  const auto& cl2Reg = this->getCluster(j).getRegion();\n"
		   "  Composition lo1 = cl1Reg.getOrigin();\n"
		   "  Composition lo2 = cl2Reg.getOrigin();\n"
		   "\n"
		   "  auto& subpaving = this->getSubpaving();\n"
		   "\n"
		   "  auto matchSpeciesPair =\n"
		   "    [&species](Composition lo, Species s1, Species s2) {\n"
		   "      if (lo[s1] == 0 || lo[s2] == 0) { return false; }\n"
		   "      for (auto l : species) {\n"
		   "        if (l == s1 || l == s2) { continue; }\n"
		   "        if (lo[l] != 0) { return false; }\n"
		   "      }\n"
		   "      return true;\n"
		   "    };\n"
		   "\n"
		   "  auto tryAddProdReaction = [&](const auto& comp) {\n"
		   "    if (comp == Composition::zero()) {\n"
		   "      return;\n"
		   "    }\n"
		   "    auto pId = subpaving.findTileId(comp);\n"
		   "    if (pId != subpaving.invalidIndex()) {\n"
		   "      this->addProductionReaction(tag, {i, j, pId});\n"
		   "    }\n"
		   "  };\n"
		   "  auto tryAddDissReaction = [&](const auto& comp) {\n"
		   "    if (comp == Composition::zero()) {\n"
		   "      return;\n"
		   "    }\n"
		   "    auto rId = subpaving.findTileId(comp);\n"
		   "    if (rId != subpaving.invalidIndex()) {\n"
		   "      this->addDissociationReaction(tag, {rId, i, j});\n"
		   "    }\n"
		   "  };\n"
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
		auto r1 = grp.parts.reactants[0];
		auto r2 = grp.parts.reactants[1];
		const auto& products = grp.parts.products;
		if (products.size() == 1 && products[0].species == nullptr) {
			// interaction cases
			ofs << "    {\n"
				   "    Composition comp = Composition::zero();\n";
			if (r1.species && r2.species) {
				auto r1SpRef = "Species::" + r1.label;
				auto r2SpRef = "Species::" + r2.label;
				ofs << "    if (lo1.isOnAxis(" << r1SpRef << ") &&\n"
					<< "        lo2.isOnAxis(" << r2SpRef << ")) {\n"
					<< "      comp[" << r1SpRef << "] =\n"
					<< "        lo1[" << r1SpRef << "];\n"
					<< "      comp[" << r2SpRef << "] =\n"
					<< "        lo2[" << r2SpRef << "];\n"
					<< "      tryAddProdReaction(comp);\n"
					<< "    }\n"
					<< "    if (lo2.isOnAxis(" << r1SpRef << ") &&\n"
					<< "        lo1.isOnAxis(" << r2SpRef << ")) {\n"
					<< "      comp[" << r1SpRef << "] =\n"
					<< "        lo2[" << r1SpRef << "];\n"
					<< "      comp[" << r2SpRef << "] =\n"
					<< "        lo1[" << r2SpRef << "];\n"
					<< "      tryAddProdReaction(comp);\n"
					<< "    }\n";
			}
			else if (xOr(r1.species, r2.species)) {
				std::string pfx = "Species::";
				auto [r1SpRef, r2Sp1Ref, r2Sp2Ref] = r1.species ?
					std::make_tuple(pfx + r1.label,
						pfx + r2.constituentLabels[0],
						pfx + r2.constituentLabels[1]) :
					std::make_tuple(pfx + r2.label,
						pfx + r1.constituentLabels[0],
						pfx + r1.constituentLabels[1]);
				ofs << "    if (lo1.isOnAxis(" << r1SpRef << ") &&\n"
					<< "        matchSpeciesPair(lo2,\n"
					<< "          " << r2Sp1Ref << ", " << r2Sp2Ref << ")) {\n"
					<< "      comp = Composition::zero();\n"
					<< "      comp[" << r1SpRef << "] +=\n"
					<< "        lo1[" << r1SpRef << "];\n"
					<< "      comp[" << r2Sp1Ref << "] +=\n"
					<< "        lo2[" << r2Sp1Ref << "];\n"
					<< "      comp[" << r2Sp2Ref << "] +=\n"
					<< "        lo2[" << r2Sp2Ref << "];\n"
					<< "      tryAddProdReaction(comp);\n"
					<< "    }\n"
					<< "    if (lo2.isOnAxis(" << r1SpRef << ") &&\n"
					<< "        matchSpeciesPair(lo1,\n"
					<< "          " << r2Sp1Ref << ", " << r2Sp2Ref << ")) {\n"
					<< "      comp = Composition::zero();\n"
					<< "      comp[" << r1SpRef << "] +=\n"
					<< "        lo2[" << r1SpRef << "];\n"
					<< "      comp[" << r2Sp1Ref << "] +=\n"
					<< "        lo1[" << r2Sp1Ref << "];\n"
					<< "      comp[" << r2Sp2Ref << "] +=\n"
					<< "        lo1[" << r2Sp2Ref << "];\n"
					<< "      tryAddProdReaction(comp);\n"
					<< "    }\n";
			}
			else {
				auto sp1Label = r1.constituentLabels[0];
				auto sp2Label = r1.constituentLabels[1];
				if (sp1Label != r2.constituentLabels[0] ||
					sp2Label != r2.constituentLabels[1]) {
					XOLOTL_ERROR(
						std::runtime_error, "unsupported reaction case");
				}
				auto sp1Ref = "Species::" + sp1Label;
				auto sp2Ref = "Species::" + sp2Label;
				ofs << "    if (matchSpeciesPair(lo1, sp1Ref, sp2Ref) &&\n"
					<< "        matchSpeciesPair(lo2, sp1Ref, sp2Ref)) {\n"
					<< "      comp[" << sp1Ref << "] +=\n"
					<< "        lo1[" << sp1Ref << "] +\n"
					<< "        lo2[" << sp1Ref << "];\n"
					<< "      comp[" << sp2Ref << "] +=\n"
					<< "        lo1[" << sp2Ref << "] +\n"
					<< "        lo2[" << sp2Ref << "];\n"
					<< "      tryAddProdReaction(comp);\n"
					<< "    }\n";
			}
			ofs << "    }\n";
			continue;
		}
		if (r1.label == r2.label) {
			auto spec = r1.label;
			auto specRef = "Species::" + spec;
			auto pSpRef = "Species::" + products[0].label;
			ofs << "    if (lo1.isOnAxis(" << specRef << ") &&\n"
				<< "        lo2.isOnAxis(" << specRef << ")) {\n"
				<< "      auto size = lo1[" << specRef << "] +\n"
				<< "        lo2[" << specRef << "];\n"
				<< "      Composition comp = Composition::zero();\n"
				<< "      comp[" << pSpRef << "] = size;\n"
				<< "      tryAddProdReaction(comp);\n"
				<< "    }\n";
		}
		else if (auto [f, v, i] = checkVacItrPair(r1.species, r2.species); f) {
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
				<< "      int prodSize = vSize - iSize;\n";
			if (products.empty()) {
				ofs << "      if (prodSize == 0) {\n"
					   "        this->addProductionReaction(tag, {i, j});\n"
					   "      }\n";
			}
			else if (products.size() == 1 &&
				products[0].species->type == SpeciesData::Type::vacancy) {
				ofs << "      if (prodSize > 0) {\n"
					<< "        Composition comp = Composition::zero();\n"
					<< "        comp[" << vSpecRef << "] = prodSize;\n"
					<< "        tryAddProdReaction(comp);\n"
					   "      }\n";
			}
			else if (products.size() == 1 &&
				products[0].species->type == SpeciesData::Type::interstitial) {
				ofs << "      if (prodSize < 0) {\n"
					<< "        Composition comp = Composition::zero();\n"
					<< "        comp[" << iSpecRef << "] = -prodSize;\n"
					<< "        tryAddProdReaction(comp);\n"
					   "      }\n";
			}
			else {
				XOLOTL_ERROR(std::runtime_error, "unsupported reaction case");
			}
			ofs << "    }\n";
		}
		else {
			XOLOTL_ERROR(std::runtime_error, "unsupported reaction case");
		}
	}

	for (auto&& [key, grp] : _dissociationReactionGroups) {
		auto p1 = grp.parts.products[0];
		auto p2 = grp.parts.products[1];
		const auto& reactants = grp.parts.reactants;
		if (reactants.size() == 1 && reactants[0].species == nullptr) {
			// interaction cases
			ofs << "    {\n"
				   "    Composition comp = Composition::zero();\n";
			if (p1.species && p2.species) {
				auto p1SpRef = "Species::" + p1.label;
				auto p2SpRef = "Species::" + p2.label;
				ofs << "    if (lo1.isOnAxis(" << p1SpRef << ") &&\n"
					<< "        lo2.isOnAxis(" << p2SpRef << ")) {\n"
					<< "      comp[" << p1SpRef << "] =\n"
					<< "        lo1[" << p1SpRef << "];\n"
					<< "      comp[" << p2SpRef << "] =\n"
					<< "        lo2[" << p2SpRef << "];\n"
					<< "      tryAddDissReaction(comp);\n"
					<< "    }\n"
					<< "    if (lo2.isOnAxis(" << p1SpRef << ") &&\n"
					<< "        lo1.isOnAxis(" << p2SpRef << ")) {\n"
					<< "      comp[" << p1SpRef << "] =\n"
					<< "        lo2[" << p1SpRef << "];\n"
					<< "      comp[" << p2SpRef << "] =\n"
					<< "        lo1[" << p2SpRef << "];\n"
					<< "      tryAddDissReaction(comp);\n"
					<< "    }\n";
			}
			else if (xOr(p1.species, p2.species)) {
				std::string pfx = "Species::";
				auto [p1SpRef, p2Sp1Ref, p2Sp2Ref] = p1.species ?
					std::make_tuple(pfx + p1.label,
						pfx + p2.constituentLabels[0],
						pfx + p2.constituentLabels[1]) :
					std::make_tuple(pfx + p2.label,
						pfx + p1.constituentLabels[0],
						pfx + p1.constituentLabels[1]);
				ofs << "    if (lo1.isOnAxis(" << p1SpRef << ") &&\n"
					<< "        matchSpeciesPair(lo2,\n"
					<< "          " << p2Sp1Ref << ", " << p2Sp2Ref << ")) {\n"
					<< "      comp = Composition::zero();\n"
					<< "      comp[" << p1SpRef << "] +=\n"
					<< "        lo1[" << p1SpRef << "];\n"
					<< "      comp[" << p2Sp1Ref << "] +=\n"
					<< "        lo2[" << p2Sp1Ref << "];\n"
					<< "      comp[" << p2Sp2Ref << "] +=\n"
					<< "        lo2[" << p2Sp2Ref << "];\n"
					<< "      tryAddDissReaction(comp);\n"
					<< "    }\n"
					<< "    if (lo2.isOnAxis(" << p1SpRef << ") &&\n"
					<< "        matchSpeciesPair(lo1,\n"
					<< "          " << p2Sp1Ref << ", " << p2Sp2Ref << ")) {\n"
					<< "      comp = Composition::zero();\n"
					<< "      comp[" << p1SpRef << "] +=\n"
					<< "        lo2[" << p1SpRef << "];\n"
					<< "      comp[" << p2Sp1Ref << "] +=\n"
					<< "        lo1[" << p2Sp1Ref << "];\n"
					<< "      comp[" << p2Sp2Ref << "] +=\n"
					<< "        lo1[" << p2Sp2Ref << "];\n"
					<< "      tryAddDissReaction(comp);\n"
					<< "    }\n";
			}
			else {
				auto sp1Label = p1.constituentLabels[0];
				auto sp2Label = p1.constituentLabels[1];
				if (sp1Label != p2.constituentLabels[0] ||
					sp2Label != p2.constituentLabels[1]) {
					XOLOTL_ERROR(
						std::runtime_error, "unsupported reaction case");
				}
				auto sp1Ref = "Species::" + sp1Label;
				auto sp2Ref = "Species::" + sp2Label;
				ofs << "    if (matchSpeciesPair(lo1,\n"
					<< "          " << sp1Ref << ", " << sp2Ref << ") &&\n"
					<< "        matchSpeciesPair(lo2,\n"
					<< "          " << sp1Ref << ", " << sp2Ref << ")) {\n"
					<< "      comp[" << sp1Ref << "] +=\n"
					<< "        lo1[" << sp1Ref << "] +\n"
					<< "        lo2[" << sp1Ref << "];\n"
					<< "      comp[" << sp2Ref << "] +=\n"
					<< "        lo1[" << sp2Ref << "] +\n"
					<< "        lo2[" << sp2Ref << "];\n"
					<< "      tryAddDissReaction(comp);\n"
					<< "    }\n";
			}
			ofs << "    }\n";
			continue;
		}
		if (reactants.size() == 1 && !(p1.species && p2.species)) {
			// One of the products is an interaction
			auto rSpRef = "Species::" + reactants[0].label;
			ofs << "    Composition comp = Composition::zero();\n";
			if (xOr(p1.species, p2.species)) {
				std::string pfx = "Species::";
				auto [p1SpRef, p2Sp1Ref, p2Sp2Ref] = p1.species ?
					std::make_tuple(pfx + p1.label,
						pfx + p2.constituentLabels[0],
						pfx + p2.constituentLabels[1]) :
					std::make_tuple(pfx + p2.label,
						pfx + p1.constituentLabels[0],
						pfx + p1.constituentLabels[1]);
				ofs << "    if (lo1.isOnAxis(" << p1SpRef << ") &&\n"
					<< "        matchSpeciesPair(lo2,\n"
					<< "          " << p2Sp1Ref << ", " << p2Sp2Ref << ")) {\n"
					<< "      comp = Composition::zero();\n"
					<< "      comp[" << p1SpRef << "] +=\n"
					<< "        lo1[" << p1SpRef << "];\n"
					<< "      comp[" << p2Sp1Ref << "] +=\n"
					<< "        lo2[" << p2Sp1Ref << "];\n"
					<< "      comp[" << p2Sp2Ref << "] +=\n"
					<< "        lo2[" << p2Sp2Ref << "];\n"
					<< "      tryAddDissReaction(comp);\n"
					<< "    }\n"
					<< "    if (lo2.isOnAxis(" << p1SpRef << ") &&\n"
					<< "        matchSpeciesPair(lo1,\n"
					<< "          " << p2Sp1Ref << ", " << p2Sp2Ref << ")) {\n"
					<< "      comp = Composition::zero();\n"
					<< "      comp[" << p1SpRef << "] +=\n"
					<< "        lo2[" << p1SpRef << "];\n"
					<< "      comp[" << p2Sp1Ref << "] +=\n"
					<< "        lo1[" << p2Sp1Ref << "];\n"
					<< "      comp[" << p2Sp2Ref << "] +=\n"
					<< "        lo1[" << p2Sp2Ref << "];\n"
					<< "      tryAddDissReaction(comp);\n"
					<< "    }\n";
			}
			continue;
		}
		if (reactants.size() == 1 && p1.label == p2.label) {
			auto spec = p1.label;
			auto specRef = "Species::" + spec;
			auto rSpRef = "Species::" + reactants[0].label;
			ofs << "    if (lo1.isOnAxis(" << specRef << ") &&\n"
				<< "        lo2.isOnAxis(" << specRef << ")) {\n"
				<< "      auto size = lo1[" << specRef << "] +\n"
				<< "        lo2[" << specRef << "];\n"
				<< "      Composition comp = Composition::zero();\n"
				<< "      comp[" << rSpRef << "] = size;\n"
				<< "      tryAddDissReaction(comp);\n"
				<< "    }\n";
			continue;
		}
		XOLOTL_ERROR(std::runtime_error, "unsupported reaction case");
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
		   "    this->_clusterData.gridSize, "
		   "this->_clusterData.numClusters,\n"
		   "    this->_enableReadRates, this->getProductionReactions(),\n"
		   "    this->getDissociationReactions(), "
		   "this->getSinkReactions());\n"
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
		<< "  ::getTotalAtomConcentration(ConcentrationsView "
		   "concentrations,\n"
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
		   "      const Region& clReg = "
		   "clData().getCluster(i).getRegion();\n"
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
		   "  auto globalData = std::vector<double>(localData.size(), "
		   "0.0);\n"
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
	ofs << "#include <plsm/EnumIndexed.h>\n"
		<< "#include <xolotl/factory/network/NetworkHandlerFactory.h>\n"
		<< "#include <" << _networkHandler << ".h>\n"
		<< "#include <" << _reactionNetwork << ".h>\n"
		<< "\n"
		   "namespace xolotl::core::network\n"
		   "{\n"
		   "namespace detail\n"
		   "{\n"
		   "using NHF = "
		   "::xolotl::factory::network::NetworkHandlerFactory;\n"
		   "template <typename T>\n"
		   "using RegCol = NHF::RegistrationCollection<T>;\n"
		   "auto networkHandlerRegistrations =\n"
		<< "  RegCol<" << _networkHandler << ">(\n"
		<< "    {\"" << _userMaterial << "\"});\n"
		<< "}\n"
		   "\n"
		   "auto nwGenerator = [](const options::IOptions& options) {\n"
		<< "  using NetworkType = " << _reactionNetwork << ";\n"
		<< "  using Species = NetworkType::Species;\n"
		   "  using AmountType = NetworkType::AmountType;\n"
		   "  const auto& params = options.getNetworkParameters();\n"
		   "  auto nSpecies = NetworkType::getNumberOfSpecies();\n"
		   "  auto maxSpeciesAmounts = std::vector<AmountType>(nSpecies);\n"
		   "  plsm::EnumIndexed<NetworkType::SubdivisionRatio, Species> rat;\n";
	for (auto&& sp : _speciesData) {
		switch (sp.type) {
		case SpeciesData::Type::impurity:
			ofs << "  rat[Species::" << sp.label << "] = params[0] + 1;\n";
			break;
		case SpeciesData::Type::vacancy:
			ofs << "  rat[Species::" << sp.label << "] = params[3] + 1;\n";
			break;
		case SpeciesData::Type::interstitial:
			ofs << "  rat[Species::" << sp.label << "] = params[4] + 1;\n";
			break;
		}
	}
	ofs << "  std::vector<NetworkType::SubdivisionRatio> subdivRatios{{rat}};\n"
		<< "  for (std::size_t i = 0; i < nSpecies; ++i) {\n"
		   "    maxSpeciesAmounts[i] = rat[i] - 1;\n"
		   "  }\n"
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
		   "using MHF = "
		   "::xolotl::factory::material::MaterialHandlerFactory;\n"
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
	fs::remove_all(_buildDir);
	fs::create_directories(_buildDir);
	XOLOTL_LOG_XTRA << "XNGEN: xolotl install prefix: " << xolotlPrefix;
	std::stringstream cmdss;
	cmdss << "cmake -S " << _genDir << " -B " << _buildDir
		  << " -DCMAKE_PREFIX_PATH=" << xolotlPrefix;

	auto outputFile = _buildDir / "out.txt";
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
