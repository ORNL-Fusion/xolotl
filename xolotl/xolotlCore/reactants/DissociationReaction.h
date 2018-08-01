#ifndef XCORE_DISSOCIATION_REACTION_H
#define XCORE_DISSOCIATION_REACTION_H

#include "ProductionReaction.h"

namespace xolotlCore {

/**
 * This is a public class that is used to store a reaction.
 *
 * The constant k- is stored along the clusters taking part in the
 * production for faster computation because they only change
 * when the temperature change. k is computed when setTemperature() is called.
 */
// Using a std::array as a key instead of a std::vector gives much
// better performance, presumably because its size is known at compile time
// instead of dynamically determined/managed.
class DissociationReaction: public KeyedReaction<
		std::array<IReactant::SizeType, 3 * NumSpecies>> {
public:

	/**
	 * The dissociating cluster
	 */
	IReactant& dissociating;

	/**
	 * The revert reaction corresponding to this dissociation
	 */
	ProductionReaction * reverseReaction;

	/**
	 * Default constructor, deleted to require construction from reactant pair.
	 */
	DissociationReaction() = delete;

	/**
	 * Construct a DissociationReaction.
	 *
	 * @param _diss Dissociating reactant to use.
	 * @param _r1 One of the reactants.
	 * @param _r2 The other reactant.
	 */
	DissociationReaction(IReactant& _diss, IReactant& _r1, IReactant& _r2,
			ProductionReaction* _reverse = nullptr) :
			KeyedReaction(_r1, _r2), dissociating(_diss), reverseReaction(
					_reverse) {

		// Build our descriptive key.
		// Assumes our first and second reactants are ordered by composition.
		auto const& dissComp = dissociating.getComposition();
		auto const& firstComp = first.getComposition();
		auto const& secondComp = second.getComposition();

		// Determine our descriptive key and cache it.
		// This assumes that our reactants are ordered.
		auto nextBegin = std::copy(dissComp.begin(), dissComp.end(),
				descKey.begin());
		nextBegin = std::copy(firstComp.begin(), firstComp.end(), nextBegin);
		std::copy(secondComp.begin(), secondComp.end(), nextBegin);
	}

	/**
	 * Copy constructor, deleted to ensure we are constructed with reactants.
	 */
	DissociationReaction(const Reaction& other) = delete;

};

/**
 * Output a basic description of the reaction to the given output stream.
 *
 * @param os Output stream on which to write description.
 * @param reaction Reaction to describe.
 * @return Output stream after having written description.
 */
std::ostream& operator<<(std::ostream& os,
		const DissociationReaction& reaction);

} // namespace xolotlCore

// For the Reaction KeyTypes  to be used as a key in an std::unordered_map,
// we need to define a hash function for it.
// Since std::unordered_map uses std::hash on its keys, and because
// ours is a user-defined type, we add our hash function to the std namespace.
//
// TODO can't this be moved up to be a templated struct on KeyedReaction,
// rather than the most derived class?  clang++ doesn't seem to like it.
namespace std {

template<>
struct hash<xolotlCore::DissociationReaction::KeyType> {

	size_t operator()(
			const xolotlCore::DissociationReaction::KeyType& val) const {
		// This may not be a good hash function - needs to be evaluated
		// for the compositions Xolotl uses.
		return std::accumulate(val.begin(), val.end(), 0);
	}
};
}

#endif /* XCORE_DISSOCIATION_REACTION_H */
