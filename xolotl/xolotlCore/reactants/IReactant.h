#ifndef IREACTANT_H
#define IREACTANT_H

// Includes
#include "IReactionNetwork.h"
#include "ReactantUtils.h"
#include <memory>
#include <vector>
#include <map>

namespace xolotlCore {

class IReactionNetwork;

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
class IReactant {

public:

	/**
	 * The destructor
	 */
	virtual ~IReactant() {
	}

	/**
	 * Returns a reactant created using the copy constructor
	 */
	virtual std::shared_ptr<IReactant> clone() = 0;

	/**
	 * Create a production pair associated with the given reaction.
	 * Create the connectivity.
	 *
	 * @param reaction The reaction creating this cluster.
	 */
	virtual void createProduction(
			std::shared_ptr<ProductionReaction> reaction) = 0;

	/**
	 * Create a combination associated with the given reaction.
	 * Create the connectivity.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 */
	virtual void createCombination(
			std::shared_ptr<ProductionReaction> reaction) = 0;

	/**
	 * Create a dissociation pair associated with the given reaction.
	 * Create the connectivity.
	 *
	 * @param reaction The reaction creating this cluster.
	 */
	virtual void createDissociation(
			std::shared_ptr<DissociationReaction> reaction) = 0;

	/**
	 * Create an emission pair associated with the given reaction.
	 * Create the connectivity.
	 *
	 * @param reaction The reaction where this cluster emits.
	 */
	virtual void createEmission(
			std::shared_ptr<DissociationReaction> reaction) = 0;

	/**
	 * Add the reactions to the network lists.
	 */
	virtual void optimizeReactions() = 0;

	/**
	 * This operation returns the current concentration.
	 *
	 * @param distA The first distance for super clusters
	 * @param distB The second distance for super clusters
	 * @return The concentration of this reactant
	 */
	virtual double getConcentration(double distA = 0.0,
			double distB = 0.0) const = 0;

	/**
	 * This operation sets the concentration of the reactant to the
	 * specified amount.
	 *
	 * @param conc The new concentation
	 */
	virtual void setConcentration(double conc) = 0;

	/**
	 * This operation returns the total flux of this reactant in the
	 * current network.
	 *
	 * @return The total change in flux for this reactant due to all
	 * reactions
	 */
	virtual double getTotalFlux() = 0;

	/**
	 * This operation sets the collection of other reactants that make up
	 * the reaction network in which this reactant exists.
	 *
	 * @param network The reaction network of which this reactant is a part
	 */
	virtual void setReactionNetwork(
			std::shared_ptr<IReactionNetwork> reactionNetwork) = 0;

	/**
	 * Release the reaction network object.
	 *
	 * This should only be done when the reaction network is no longer needed
	 * by the program, and is done to break dependence cycles that would
	 * otherwise keep the network and reactant objects from being destroyed.
	 */
	virtual void releaseReactionNetwork() = 0;

	/**
	 * This operation signifies that the reactant with reactant Id should be
	 * listed as connected with this reactant through forward reactions.
	 *
	 * @param id The integer id of the reactant that is connected
	 * to this reactant
	 */
	virtual void setReactionConnectivity(int id) = 0;

	/**
	 * This operation signifies that the reactant with reactant Id should be
	 * listed as connected with this reactant through forward reactions.
	 *
	 * @param id The integer id of the reactant that is connected
	 * to this reactant
	 */
	virtual void setDissociationConnectivity(int id) = 0;

	/**
	 * This operation reset the connectivity sets based on the information
	 * in the effective production and dissociation vectors.
	 */
	virtual void resetConnectivities() = 0;

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
	virtual std::vector<int> getConnectivity() const = 0;

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
	virtual std::vector<double> getPartialDerivatives() const = 0;

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
	virtual void getPartialDerivatives(
			std::vector<double> & partials) const = 0;

	/**
	 * This operation returns the name of the reactant.
	 *
	 * @return The name
	 */
	virtual const std::string getName() const = 0;

	/**
	 * This operation returns the reactant's type. It is up to subclasses to
	 * define exactly what the allowed types may be.
	 *
	 * @return The type of this reactant as a string
	 */
	virtual std::string getType() const = 0;

	/**
	 * This operation returns the composition of this reactant. This map is empty
	 * when returned by the base class.
	 *
	 * @return The composition returned as a map with keys naming distinct
	 * elements and values indicating the amount of the element present.
	 */
	virtual const std::map<std::string, int> & getComposition() const = 0;

	/**
	 * Get a string containing the canonical representation of the
	 * composition of this reactant.  The string is not intended to
	 * be human-readable, but rather is useful for keys in reactant maps
	 * and for composition match tests (as opposed to comparisons of
	 * the composition maps themselves).
	 * TODO is this the same information as our name?
	 *
	 * @return A string containing the canonical representation of our
	 * composition.
	 */
	virtual std::string getCompositionString() const = 0;

	/**
	 * This operation sets the id of the reactant, The id is zero by default
	 * and clients, most likely the ReactionNetwork, are expected to set the
	 * id as needed.
	 *
	 * @param nId The new id for this reactant
	 */
	virtual void setId(int nId) = 0;

	/**
	 * This operation returns the id for this reactant.
	 *
	 * @return The id
	 */
	virtual int getId() const = 0;

	/**
	 * This operation sets the id of the moment of the reactant.
	 *
	 * @param nId The new id for this moment
	 */
	virtual void setMomentId(int nId) = 0;

	/**
	 * This operation returns the id for this reactant moment.
	 *
	 * @return The id
	 */
	virtual int getMomentId() const = 0;

	/**
	 * This operation sets the id of the helium moment of the reactant.
	 *
	 * @param nId The new id for this moment
	 */
	virtual void setHeMomentId(int nId) = 0;

	/**
	 * This operation returns the id for this reactant helium moment.
	 *
	 * @return The id
	 */
	virtual int getHeMomentId() const = 0;

	/**
	 * This operation sets the id of the vacancy moment of the reactant.
	 *
	 * @param nId The new id for this moment
	 */
	virtual void setVMomentId(int nId) = 0;

	/**
	 * This operation returns the id for this reactant vacancy moment.
	 *
	 * @return The id
	 */
	virtual int getVMomentId() const = 0;

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
	 * @param temp The new reactant temperature
	 */
	virtual void setTemperature(double temp) = 0;

	/**
	 * This operation returns the temperature at which the reactant currently exists.
	 *
	 * @return The temperature.
	 */
	virtual double getTemperature() const = 0;

	/**
	 * This operation returns the total size of the reactant.
	 *
	 * @return The total size of this reactant including the contributions
	 * from all species types
	 */
	virtual int getSize() const = 0;

	/**
	 * This operation retrieves the formation energy for this reactant.
	 *
	 * @return The value of the formation energy
	 */
	virtual double getFormationEnergy() const = 0;

	/**
	 * This operation sets the formation energy for this reactant.
	 *
	 * @param energy The formation energy
	 */
	virtual void setFormationEnergy(double energy) = 0;

	/**
	 * This operation retrieves the diffusion factor, D_0, that is used to
	 * calculate the diffusion coefficient for this reactant.
	 *
	 * @return The diffusion factor of this reactant
	 */
	virtual double getDiffusionFactor() const = 0;

	/**
	 * This operation sets the diffusion factor, D_0, that is used to calculate
	 * the diffusion coefficient for this reactant.
	 *
	 * @param factor The diffusion factor
	 */
	virtual void setDiffusionFactor(const double factor) = 0;

	/**
	 * This operation returns the diffusion coefficient for this reactant and is
	 * calculated from the diffusion factor.
	 *
	 * @return The diffusion coefficient
	 */
	virtual double getDiffusionCoefficient() const = 0;

	/**
	 * This operation sets the migration energy for this reactant.
	 *
	 * @param energy The migration energy
	 */
	virtual void setMigrationEnergy(const double energy) = 0;

	/**
	 * This operation retrieves the migration energy for this reactant.
	 *
	 * @return the migration energy
	 */
	virtual double getMigrationEnergy() const = 0;

	/**
	 * This operation returns the reaction radius for the
	 * particular reactant.
	 *
	 * @return The reaction radius
	 */
	virtual double getReactionRadius() const = 0;

	/**
	 * This operation returns the sum of combination rate and emission rate
	 * (where this reactant is on the left side of the reaction) for this
	 * particular reactant.
	 * This is used to computed the desorption rate in the
	 * modified trap-mutation handler.
	 *
	 * @return The rate
	 */
	virtual double getLeftSideRate() const = 0;

	/**
	 * This operation returns true if the cluster is a mixed-species or compound
	 * cluster and false if it is a single species cluster.
	 */
	virtual bool isMixed() const = 0;

};

} // end namespace xolotlCore

#endif
