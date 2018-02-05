#ifndef XCORE_PENDING_PRODUCTION_REACTION_INFO_H
#define XCORE_PENDING_PRODUCTION_REACTION_INFO_H

namespace xolotlCore {

class IReactant;

/**
 * Information about a production reaction that needs to be created
 * for an unspecified pair of reactants.
 * Used to support batch creation of production and dissociation reactions
 * for a reactant interacting with a super cluster.
 */
struct PendingProductionReactionInfo {

	IReactant& product;
	int numHe;
	int numV;
	int i;
	int j;

	PendingProductionReactionInfo(IReactant& _product, int _numHe, int _numV,
			int _i, int _j) :
			product(_product), numHe(_numHe), numV(_numV), i(_i), j(_j) {
	}

	/**
	 * Default and copy ctor, disallowed to detect potential use.
	 */
	PendingProductionReactionInfo() = delete;

	/**
	 * Copy ctor, using default implementation.
	 */
	PendingProductionReactionInfo(const PendingProductionReactionInfo& other) = default;

	/**
	 * Move ctor, using default implementation.  Needed if
	 * PendingProductionReactionInfos are stored in vectors.
	 */
	PendingProductionReactionInfo(PendingProductionReactionInfo&& other) = default;
};

} // namespace xolotlCore

#endif /* XCORE_PENDING_PRODUCTION_REACTION_INFO_H */
