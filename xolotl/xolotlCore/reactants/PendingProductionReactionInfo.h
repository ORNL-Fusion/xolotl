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
	int a[4] = { };
	int b[4] = { };

	PendingProductionReactionInfo(IReactant& _product, int _numA[4] = { },
			int _numB[4] = { }) :
			product(_product) {
		a[0] = _numA[0], a[1] = _numA[1], a[2] = _numA[2], a[3] = _numA[3];
		b[0] = _numB[0], b[1] = _numB[1], b[2] = _numB[2], b[3] = _numB[3];
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
