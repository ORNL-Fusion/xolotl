#ifndef REACTANT_H
#define REACTANT_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include <map>
#include "ReactionNetwork.h"

namespace xolotlCore {

/**
 * A Reactant is a general reacting body in a reaction network. It represents
 * any body whose population can change with time due to reactions of any type.
 *
 * Reactants inherently know the other Reactants with which they interact. They
 * declare their interactions with other Reactants in the network after it is
 * set (setReactionNetwork) via the getConnectivity() operation. "Connectivity"
 * indicates whether two Reacants interact, via any mechanism, in an abstract
 * sense (as if they were nodes connected by an edge on a network graph).
 *
 * This is an abstract base class that only provides direct support for
 * manipulate the concentration, etc. It should be subclassed to add
 * functionality for calculate fluxes and computing connectivity.
 */
class Reactant {

protected:

	/** The total concentration of this Reactant.
	 */
	double concentration;

	/** The name of this Reactant.
	 */
	std::string name;

	/** An integer identification number for this reactant.
	 */
	int id;

	/** The reaction network that includes this reactant.
	 */
	std::shared_ptr<ReactionNetwork> network;

	/**
	 * The map that contains the composition of this cluster
	 */
	std::map<std::string,int> compositionMap;

public:

	/** The constructor.
	 */
	Reactant();

	/**
	 * The copy constructor. All reactants MUST be deep copied.
	 * @param other The reactant to copy
	 */
	Reactant(const Reactant &other);

	/** The destructor
	 */
	virtual ~Reactant();

	/**
	 * This operation returns a Reactant that is created using the copy
	 * constructor. If this Reactant is actually a subclass of Reactant, the
	 * clone will be of the same type and therefore carry all of the members
	 * and virtual functions of the subclass in addition to those of the
	 * Reactant. This type of copy is not only handy but, in fact, quite
	 * necessary in those cases where a Reactant must be copied but its exact
	 * subclass is unknown and there is no way to make a reasonable assumption
	 * about it.
	 * @return A copy of this reactant.
	 */
	virtual std::shared_ptr<Reactant> clone();

	/**
	 * An alternative constructor that can be used to create a reactant
	 * with an initial concentration.
	 *
	 * @param conc The initial concentration
	 */
	Reactant(double conc);

	/**
	 * This operation returns the current concentration.
	 *
	 * @return The concentration of this reactant
	 */
	double getConcentration() const;

	/**
	 * This operation increases the concentration of the reactant by the
	 * specified amount.
	 *
	 * @param deltaConc the change in concentration
	 */
	void increaseConcentration(double deltaConc);

	/**
	 * This operation decreases the concentration of the reactant by the
	 * specified amount.
	 *
	 * @param deltaConc the change in concentration
	 */
	void decreaseConcentration(double deltaConc);

	/**
	 * This operation sets the concentration of the reactant to the
	 * specified amount.
	 *
	 * @param conc The new concentation
	 */
	void setConcentration(double conc);

	/**
	 * This operation sets the concentration of the reactant to zero.
	 */
	void zero();

	/**
	 * This operation returns the total flux of this reactant in the
	 * current network.
	 *
	 * @param temperature The temperature at which to calculate the Diffusion Coefficient
	 * @return The total change in flux for this reactant due to all
	 * reactions
	 */
	virtual double getTotalFlux(const double temperature);

	/**
	 * This operation sets the collection of other reactants that make up
	 * the reaction network in which this reactant exists.
	 *
	 * @param network The reaction network of which this reactant is a part
	 */
	virtual void setReactionNetwork(
			std::shared_ptr<ReactionNetwork> reactionNetwork);

	/**
	 * This operation returns a list that represents the connectivity
	 * between this Reactant and other Reactants in the network.
	 * "Connectivity" indicates whether two Reactants interact, via any
	 * mechanism, in an abstract sense (as if they were nodes connected by
	 * an edge on a network graph).
	 *
	 * @return An array of ones and zeros that indicate whether or not this
	 * Reactant interacts via any mechanism with another Reactant. A "1" at
	 * the i-th entry in this array indicates that the Reactant interacts
	 * with the i-th Reactant in the ReactionNetwork and a "0" indicates
	 * that it does not.
	 */
	virtual std::vector<int> getConnectivity() const;

	/**
	 * This operation returns the list of partial derivatives of this Reactant
	 * with respect to all other reactants in the network. The combined lists
	 * of partial derivatives from all of the reactants in the network can be
	 * used to form, for example, a Jacobian.
	 *
	 * @param the temperature at which the reactions are occurring
	 * @return The partial derivatives for this reactant where index zero
	 * corresponds to the first reactant in the list returned by the
	 * ReactionNetwork::getAll() operation.
	 */
	virtual std::vector<double> getPartialDerivatives(double temperature) const;

	/**
	 * This operation writes the contents of the reactant to a string. This
	 * operation should be overridden by subclasses.
	 *
	 * @return A serialized version of this reactant as a string.
	 */
	virtual const std::string toString();

	/**
	 * This operation returns the name of the reactant.
	 * @return the name
	 */
	const std::string getName() const;

	/**
	 * This operation returns the compositon of this reactant. This map is empty
	 * when returned by the base class.
	 * @return The composition returned as a map with keys naming distinct
	 * elements and values indicating the amount of the element present.
	 */
	virtual const std::map<std::string, int> & getComposition() const;
};

} // end namespace xolotlCore

#endif
