#ifndef XOLOTL_CORE_REACTANTTYPE_H_
#define XOLOTL_CORE_REACTANTTYPE_H_

#include <string>

#include "Species.h"

namespace xolotlCore {

// TODO Don't like having to define ReactantTypes for all types of
// ReactionNetworks here.  Would prefer that the ReactantTypes are
// defined in the the derived ReactionNetwork subdirectories,
// so that ReactionNetworks can define those for only those that
// it knows about.
//
// C++ doesn't support derivation of enum classes, such that we
// could define a base set of ReactantTypes and then extend it
// by the specific types of reaction networks.
//
// Since Reactant keeps its own type, and the other parts of the code
// also use these reactant types without distinguishing whether they
// have a PSIClusterReactionNetwork or NEClusterReactionNetwork,
// and because support for both types of reaction entworks is always built,
// for now we're stuck defining all potential reactant types in this shared
// enum.
//
// We could finesse the problem using preprocessor tricks, e.g., by
// including a file from the psiclusters and neclusters directories
// in the appropriate place in the enum definition.
//
// Note: do not assume that the values of these enums correspond to
// those in the Species enum.
//
enum class ReactantType {
	// Common Reactant types.
	Invalid = -1,
	V,
	I,
	He,
	D,
	T,
	HeV,
	HeI,
	Xe,
	XeV,
	XeI,
	Void,
	Frank,
	Faulted,
	Perfect,

	// Reactant types for PSI networks
	PSIMixed,
	PSISuper,

	// Reactant types for NE networks
	NESuper,

	// Reactant types for iron networks
	FeSuper,

	// Reactant types for alloy networks
	VoidSuper,
	FrankSuper,
	FaultedSuper
};

std::string toString(ReactantType rtype);

Species toSpecies(ReactantType rtype);

/**
 * Obtain the ReactantType that corresponds to a given Species.
 *
 * @param s The Species of interest.
 * @return The ReactantType taht corresponds to Species s.
 */
ReactantType toReactantType(Species s);

} // namespace xolotlCore

// For a ReactantType to be used as a key in an std::unordered_map,
// we need a hash function for it.
// Some compilers seem to automatically handle this, some need us
// to do it explicitly.
// Since std::unordered_map uses std::hash on its keys, and because
// ours is a user-defined type, we add our hash function to the std namespace.
namespace std {

template<>
struct hash<xolotlCore::ReactantType> {

	size_t operator()(const xolotlCore::ReactantType rtype) const {
		return static_cast<size_t>(rtype);
	}
};
} // namespace std

#endif /* XOLOTL_CORE_REACTANTTYPE_H_ */
