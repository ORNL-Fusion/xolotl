#include <string>
#include <unordered_map>
#include "Species.h"

namespace xolotlCore {

std::string toString(Species s) {

	static std::unordered_map<Species, std::string> smap { { Species::Invalid,
			"Invalid_species" }, { Species::V, "V" }, { Species::I, "I" }, {
			Species::He, "He" }, { Species::D, "D" }, { Species::T, "T" }, {
			Species::Xe, "Xe" }, };

	auto iter = smap.find(s);
	return (iter != smap.end()) ? iter->second : "[unrecognized species]";
}

} // namespace xolotlCore

