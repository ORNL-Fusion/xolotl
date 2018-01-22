#ifndef REACTANTUTILS_H
#define REACTANTUTILS_H

#include <string>
#include "IReactant.h"

namespace xolotlCore {

class IReactant;

/**
 * This is a public class that is used to store a reaction.
 *
 * The constant k is stored here and other information is implemented by
 * the daughter classes. k is computed when setTemperature() is called.
 */
class Reaction {
public:
    /**
     * Type of a canonical key describing this ProductionReaction that 
     * can be used to compare it to other ProductionReactions.
     */
    typedef std::string KeyType;

protected:
    /**
     * A descriptive key in canonical form that can be used for
     * fast compares against that of other ProductionReactions.
     */
    KeyType descKey;

public:

	/**
	 * The rate constant
	 */
	double kConstant;

	//! The constructor
	Reaction() :
			kConstant(0.0) {
	}


    /**
     * Find the canonical key describing this ProductionReaction.
     */
    KeyType descriptiveKey() const  { return descKey; }
};

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
class ProductionReaction: public Reaction {
public:
	/**
	 * The first cluster in the pair.
     * Its composition string is "<=" that of the second cluster.
	 */
	IReactant * first;

	/**
	 * The second cluster in the pair
     * Its composition string is ">=" that of the second cluster.
	 */
	IReactant * second;

	//! The constructor
	ProductionReaction(IReactant * _reactant1, IReactant * _reactant2);
};

/**
 * This is a public class that is used to store a reaction.
 *
 * The constant k- is stored along the clusters taking part in the
 * production for faster computation because they only change
 * when the temperature change. k is computed when setTemperature() is called.
 */
class DissociationReaction: public Reaction {
public:

	/**
	 * The dissociating cluster
	 */
	IReactant * dissociating;

	/**
	 * The first cluster in the pair
	 */
	IReactant * first;

	/**
	 * The second cluster in the pair
	 */
	IReactant * second;

	/**
	 * The revert reaction corresponding to this dissociation
	 */
	ProductionReaction * reverseReaction;

	//! The constructor
	DissociationReaction(IReactant * dissociatingPtr, IReactant * firstPtr,
			IReactant * secondPtr);
};

}

#endif /* REACTANTUTILS_H */
