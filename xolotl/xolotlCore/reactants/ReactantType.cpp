#include <string>
#include <unordered_map>
#include "ReactantType.h"

namespace xolotlCore {

std::string toString(ReactantType rtype) {

	static std::unordered_map<ReactantType, std::string> smap { {
			ReactantType::Invalid, "Invalid_reactant_type" }, { ReactantType::V,
			"V" }, { ReactantType::I, "I" }, { ReactantType::He, "He" }, {
			ReactantType::D, "D" }, { ReactantType::T, "T" }, {
			ReactantType::HeV, "HeV" }, { ReactantType::HeI, "HeI" }, {
			ReactantType::Xe, "Xe" }, { ReactantType::XeV, "XeV" }, {
			ReactantType::XeI, "XeI" }, { ReactantType::PSIMixed, "PSIMixed" },
			{ ReactantType::PSISuper, "PSISuper" }, { ReactantType::NESuper,
					"NESuper" }, { ReactantType::FeSuper, "FeSuper" } };

	auto iter = smap.find(rtype);
	return (iter != smap.end()) ? iter->second : "[unrecognized reactant type]";
}

Species toSpecies(ReactantType rtype) {
	// Map of all single species reactant types, plus Invalid.
	static std::unordered_map<ReactantType, Species> smap { {
			ReactantType::Invalid, Species::Invalid }, { ReactantType::V,
			Species::V }, { ReactantType::I, Species::I }, { ReactantType::He,
			Species::He }, { ReactantType::D, Species::D }, { ReactantType::T,
			Species::T }, { ReactantType::Xe, Species::Xe } };

	auto iter = smap.find(rtype);
	return (iter != smap.end()) ? iter->second : Species::Invalid;
}

ReactantType toReactantType(Species s) {
	// Map of all single species reactant types, plus Invalid.
	static std::unordered_map<Species, ReactantType> smap { { Species::Invalid,
			ReactantType::Invalid }, { Species::V, ReactantType::V }, {
			Species::I, ReactantType::I }, { Species::He, ReactantType::He }, {
			Species::D, ReactantType::D }, { Species::T, ReactantType::T }, {
			Species::Xe, ReactantType::Xe } };

	auto iter = smap.find(s);
	return (iter != smap.end()) ? iter->second : ReactantType::Invalid;
}

} // namespace xolotlCore

