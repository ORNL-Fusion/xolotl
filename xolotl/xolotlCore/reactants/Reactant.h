#ifndef REACTANT_H
#define REACTANT_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include <map>
#include "ReactionNetwork.h"
#include <iostream>

namespace xolotlPerf {
	class IHandlerRegistry;
	class IEventCounter;
}

namespace xolotlCore {

/**
 * A reactant is a general reacting body in a reaction network. It represents
 * any body whose population can change with time due to reactions of any type.
 *
 * Reactants inherently know the other reactants with which they interact. They
 * declare their interactions with other reactants in the network after it is
 * set (setReactionNetwork) via the getConnectivity() operation. "Connectivity"
 * indicates whether two Reacants interact, via any mechanism, in an abstract
 * sense (as if they were nodes connected by an edge on a network graph).
 *
 * This is an abstract base class that only provides direct support for
 * manipulating the concentration, etc. It should be subclassed to add
 * functionality for calculate fluxes and computing connectivity.
 */
class Reactant {

protected:

	/**
	 * The total concentration of this reactant.
	 */
	double concentration;

	/**
	 * The name of this reactant.
	 */
	std::string name;

	/**
	 * The type name of the reactant.
	 */
	std::string typeName;

	/**
	 * An integer identification number for this reactant.
	 */
	int id;

	/**
	 * The temperature at which the cluster currently exists. The diffusion
	 * coefficient is recomputed each time the temperature is changed.
	 */
	double temperature;

	/**
	 * The reaction network that includes this reactant.
	 */
	std::shared_ptr<ReactionNetwork> network;

	/**
	 * The map that contains the composition of this cluster.
	 */
	std::map<std::string, int> compositionMap;

	/**
	 * The performance handler registry that will be used with
	 * this class.
	 */
	std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry;

	/**
	 * The constructor.
	 */
	Reactant();

public:

	/**
	 * The constructor.
	 *
	 * @param registry The performance handler registry to use
	 */
	Reactant(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * An alternative constructor that can be used to create a reactant
	 * with an initial concentration.
	 *
	 * @param conc The initial concentration
	 * @param registry The performance handler registry to use
	 */
	Reactant(double conc,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * The copy constructor. All reactants MUST be deep copied.
	 *
	 * @param other The reactant to copy
	 */
	Reactant(const Reactant &other);

	/**
	 * The destructor
	 */
	virtual ~Reactant() {}

	/**
	 * This operation returns a reactant that is created using the copy
	 * constructor. If this reactant is actually a subclass of reactant, the
	 * clone will be of the same type and therefore carry all of the members
	 * and virtual functions of the subclass in addition to those of the
	 * reactant. This type of copy is not only handy but, in fact, quite
	 * necessary in those cases where a reactant must be copied but its exact
	 * subclass is unknown and there is no way to make a reasonable assumption
	 * about it.
	 *
	 * @return A copy of this reactant.
	 */
	virtual std::shared_ptr<Reactant> clone();

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
	 * @return The total change in flux for this reactant due to all
	 * reactions
	 */
	virtual double getTotalFlux();

	/**
	 * This operation sets the collection of other reactants that make up
	 * the reaction network in which this reactant exists.
	 *
	 * @param network The reaction network of which this reactant is a part
	 */
	virtual void setReactionNetwork(
			std::shared_ptr<ReactionNetwork> reactionNetwork);

	/**
	 * Release the reaction network object.
	 *
	 * This should only be done when the reaction network is no longer needed
	 * by the program, and is done to break dependence cycles that would
	 * otherwise keep the network and reactant objects from being destroyed.
	 */
	virtual void releaseReactionNetwork() {network.reset();}

	/**
	 * This operation returns a list that represents the connectivity
	 * between this reactant and other reactants in the network.
	 * "Connectivity" indicates whether two reactants interact, via any
	 * mechanism, in an abstract sense (as if they were nodes connected by
	 * an edge on a network graph).
	 *
	 * @return An array of ones and zeros that indicate whether or not this
	 * reactant interacts via any mechanism with another reactant. A "1" at
	 * the i-th entry in this array indicates that the reactant interacts
	 * with the i-th reactant in the ReactionNetwork and a "0" indicates
	 * that it does not.
	 */
	virtual std::vector<int> getConnectivity() const;

	/**
	 * This operation returns the list of partial derivatives of this reactant
	 * with respect to all other reactants in the network. The combined lists
	 * of partial derivatives from all of the reactants in the network can be
	 * used to form, for example, a Jacobian.
	 *
	 * @return the partial derivatives for this reactant where index zero
	 * corresponds to the first reactant in the list returned by the
	 * ReactionNetwork::getAll() operation.
	 */
	virtual std::vector<double> getPartialDerivatives() const;

	/**
	 * This operation works as getPartialDerivatives above, but instead of
	 * returning a vector that it creates it fills a vector that is passed to
	 * it by the caller. This allows the caller to optimize the amount of
	 * memory allocations to just one if they are accessing the partial
	 * derivatives many times.
	 *
	 * The base class (Reactant) implementation does nothing.
	 *
	 * @param partials The vector that should be filled with the partial derivatives
	 * for this reactant where index zero corresponds to the first reactant in
	 * the list returned by the ReactionNetwork::getAll() operation. The size of
	 * the vector should be equal to ReactionNetwork::size().
	 */
	virtual void getPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation returns the name of the reactant.
	 *
	 * @return The name
	 */
	const std::string getName() const;

	/**
	 * This operation returns the reactant's type. It is up to subclasses to
	 * define exactly what the allowed types may be.
	 *
	 * @return The type of this reactant as a string
	 */
	std::string getType() const;

	/**
	 * This operation returns the composition of this reactant. This map is empty
	 * when returned by the base class.
	 *
	 * @return The composition returned as a map with keys naming distinct
	 * elements and values indicating the amount of the element present.
	 */
	virtual const std::map<std::string, int> & getComposition() const;

	/**
	 * This operation sets the id of the reactant, The id is zero by default
	 * and clients, most likely the ReactionNetwork, are expected to set the
	 * id as needed.
	 *
	 * @param nId The new id for this reactant
	 */
	void setId(int nId) {id = nId;}

	/**
	 * This operation returns the id for this reactant.
	 *
	 * @return The id
	 */
	int getId() const {return id;}

	/**
	 * This operation sets the temperature at which the reactant currently
	 * exists. Temperature-dependent quantities are recomputed when this
	 * operation is called, so the temperature should always be set first.
	 *
	 * The simplest way to set the temperature for all reactants is to call the
	 * ReactionNetwork::setTemperature() operation.
	 *
	 * The base class implementation only stores the temperature value.
	 * Subclasses are responsible for implementing their own update
	 * calculations and for calling setTemperature() in their copy constructors.
	 *
	 * @param temp The new cluster temperature
	 */
	virtual void setTemperature(double temp) {temperature = temp;}

	/**
	 * This operation returns the temperature at which the reactant currently exists.
	 *
	 * @return The temperature.
	 */
	double getTemperature() const {return temperature;}

	/**
	 * Function to overload the streaming operator in order to output reactant
	 * information easily.
	 *
	 * @param out The output stream
	 * @param reactant The reactant
	 * @return The output stream that will print desired reactant information
	 */
	friend std::ostream& operator<<(std::ostream& out,
			const Reactant& reactant);
};

} // end namespace xolotlCore

#endif
