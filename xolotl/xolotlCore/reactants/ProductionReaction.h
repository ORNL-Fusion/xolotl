#ifndef XCORE_PRODUCTION_REACTION_H
#define XCORE_PRODUCTION_REACTION_H

#include <array>
#include "Reaction.h"

namespace xolotlCore {

/**
 * This is a public class that is used to store a reaction.
 *
 * The constant k+ is stored along the clusters taking part in the
 * production for faster computation because they only change
 * when the temperature change.
 * We don't keep the product cluster because it is not needed for k calculations.
 *
 * k is computed when setTemperature() is called.
 */
// Using a std::array as a key instead of a std::vector gives much
// better performance, presumably because its size is known at compile time
// instead of dynamically determined/managed.
class ProductionReaction: public KeyedReaction<
		std::array<IReactant::SizeType,
				2 * (static_cast<int>(Species::last) + 1)> > {
public:
	/**
	 * Information about a production reaction that needs to be created
	 * for an unspecified pair of reactants.
	 * Used to support batch creation of production and dissociation reactions
	 * for a reactant interacting with a super cluster.
	 */
	struct PendingInfo {

		IReactant& product;
		int numHe;
		int numV;
		int i;
		int j;

		PendingInfo(IReactant& _product, int _numHe, int _numV, int _i, int _j) :
				product(_product), numHe(_numHe), numV(_numV), i(_i), j(_j) {
		}

		/**
		 * Default and copy ctor, disallowed to detect potential use.
		 */
		PendingInfo() = delete;
		PendingInfo(const PendingInfo& other) = delete;

		/**
		 * Move ctor, using default implementation.  Needed if
		 * PendingInfos are stored in vectors.
		 */
		PendingInfo(PendingInfo&& other) = default;
	};

	/**
	 * Default constructor, deleted to require construction from reactant pair.
	 */
	ProductionReaction() = delete;

	/**
	 * Copy constructor, deleted to ensure we are constructed with reactants.
	 */
	ProductionReaction(const Reaction& other) = delete;

	/**
	 * Construct a ProductionReaction.
	 *
	 * @param _r1 One of the reactants.
	 * @param _r2 The other reactant.
	 */
	ProductionReaction(IReactant& _r1, IReactant& _r2) :
			KeyedReaction(_r1, _r2) {

		// Build our decriptive key.
		// Assumes that our reactants are ordered by composition.
		auto const& firstComp = first.getComposition();
		auto const& secondComp = second.getComposition();

		auto nextBegin = std::copy(firstComp.begin(), firstComp.end(),
				descKey.begin());
		std::copy(secondComp.begin(), secondComp.end(), nextBegin);
	}
};

/**
 * Output a basic description of the reaction to the given output stream.
 *
 * @param os Output stream on which to write description.
 * @param reaction Reaction to describe.
 * @return Output stream after having written description.
 */
std::ostream& operator<<(std::ostream& os, const ProductionReaction& reaction);

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
struct hash<xolotlCore::ProductionReaction::KeyType> {

	size_t operator()(
			const xolotlCore::ProductionReaction::KeyType& val) const {
		// This may not be a good hash function - needs to be evaluated
		// for the compositions Xolotl uses.
		return std::accumulate(val.begin(), val.end(), 0);
	}
};
}

#endif /* XCORE_PRODUCTION_REACTION_H */
