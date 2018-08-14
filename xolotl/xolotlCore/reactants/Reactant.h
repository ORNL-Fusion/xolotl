#ifndef REACTANT_H
#define REACTANT_H

// Includes
#include <math.h>
#include <sstream>
#include <set>
#include "IReactant.h"
#include "IReactionNetwork.h"
#include "ProductionReaction.h"
#include "DissociationReaction.h"

namespace xolotlPerf {
class IHandlerRegistry;
class IEventCounter;
}

// We use std::unordered_map for quick lookup of info about
// reactions we participate in.
// The C++ standard library defines a std::hash for keys
// that are a single pointer, but not for pairs of pointers,
// so we define our own here.  To improve readability,
// we define a concise name for type of a pair of IReactant pointers
// that we use as keys.
// TODO should this be moved "upward," e.g., into IReactant.h?
namespace xolotlCore {
using ReactantAddrPair = std::pair<IReactant*, IReactant*>;
} // namespace xolotlCore

namespace std {

template<>
struct hash<xolotlCore::ReactantAddrPair> {
	size_t operator()(const xolotlCore::ReactantAddrPair& pr) const {
		// Idea for implementation taken from
		// https://www.sultanik.com/blog/HashingPointers.
		auto sum = reinterpret_cast<uintptr_t>(pr.first)
				+ reinterpret_cast<uintptr_t>(pr.second);
		// Ensure result will fit in size_t
#if SIZE_MAX < UINTPTR_MAX
		sum %= SIZE_MAX;
#endif // SIZE_MAX < UINTPTR_MAX
		return sum;
	}
};

} // namespace std

namespace xolotlCore {

/**
 * A reactant is a general reacting body in a reaction network. It represents
 * any body whose population can change with time due to reactions of any type.
 *
 * Reactants inherently know the other reactants with which they interact. They
 * declare their interactions with other reactants in the network after it is
 * set (updateFromNetwork) via the getConnectivity() operation. "Connectivity"
 * indicates whether two Reacants interact, via any mechanism, in an abstract
 * sense (as if they were nodes connected by an edge on a network graph).
 *
 * This is an abstract base class that only provides direct support for
 * manipulating the concentration, etc. It should be subclassed to add
 * functionality for calculate fluxes and computing connectivity.
 */
class Reactant: public IReactant {

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
	ReactantType type;

	/**
	 * An integer identification number for this reactant.
	 */
	int id;

	/**
	 * An integer identification number for the first moments.
	 */
	int momId[4] = { };

	/**
	 * The temperature at which the cluster currently exists. The diffusion
	 * coefficient is recomputed each time the temperature is changed.
	 */
	std::vector<double> temperature;

	/**
	 * The reaction network that includes this reactant.
	 */
	IReactionNetwork& network;

	/**
	 * The map that contains the composition of this cluster.
	 */
	IReactant::Composition composition;

	/**
	 * The performance handler registry that will be used with
	 * this class.
	 */
	std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry;

	/**
	 * The total size of this cluster including the contributions from all
	 * species.
	 */
	IReactant::SizeType size;

	/**
	 * The diffusion factor, D_0, that is used to calculate the diffusion
	 * coefficient for this cluster. The default value is 0 (does not diffuse).
	 */
	double diffusionFactor;

	/**
	 * The diffusion coefficient computed from the diffusion factor using an
	 * Arrhenius rate equation. It is re-computed every time the temperature is
	 * updated.
	 */
	std::vector<double> diffusionCoefficient;

	/**
	 * The formation energy of this cluster. It will be used to compute the
	 * binding energies appearing in the dissociation constant calculation.
	 */
	double formationEnergy;

	/**
	 * The migration energy for this cluster.
	 */
	double migrationEnergy;

	/**
	 * The reaction radius of this cluster
	 */
	double reactionRadius;

	/**
	 * The row of the reaction connectivity matrix corresponding to
	 * this Reactant stored as a set.
	 *
	 * If a cluster is involved in a reaction with this Reactant,
	 * the cluster id is an element of this set.
	 */
	std::set<int> reactionConnectivitySet;

	/**
	 * The row of the dissociation connectivity matrix corresponding to
	 * this Reactant stored as a set.
	 *
	 * If this Reactant can dissociate into a particular cluster,
	 * the cluster id is an element of this set.
	 */
	std::set<int> dissociationConnectivitySet;

	/**
	 * This operation recomputes the diffusion coefficient. It is called
	 * whenever the diffusion factor, migration energy or temperature change.
	 *
	 * @param temp the temperature
	 * @param i The position on the grid
	 */
	void recomputeDiffusionCoefficient(double temp, int i);

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	Reactant() = delete;

	/**
	 * The constructor.
	 *
	 * @param _network The network we will belong to.
	 * @param _name Our human-readable name.
	 * @param _registry The performance handler registry to use
	 */
	Reactant(IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> _registry,
			const std::string& _name = "Reactant");

	/**
	 * Copy constructor.
	 * Only used to construct dummy cluster objects of Reactant type
	 * as a copy of the Reactant part of objects of a more derived type.
	 * The more derived types initialize their base class' data, and
	 * we don't have a ctor that lets them specify all of our data,
	 * so we use this ctor to copy the Reactant data.
	 *
	 * @param other The reactant to copy
	 */
	Reactant(Reactant &other) :
			concentration(other.concentration), name(other.name), type(
					other.type), id(other.id), temperature(other.temperature), network(
					other.network), handlerRegistry(other.handlerRegistry), size(
					other.size), composition(other.composition), formationEnergy(
					other.formationEnergy), diffusionFactor(
					other.diffusionFactor), diffusionCoefficient(
					other.diffusionCoefficient), migrationEnergy(
					other.migrationEnergy), reactionRadius(
					other.reactionRadius), reactionConnectivitySet(
					other.reactionConnectivitySet), dissociationConnectivitySet(
					other.dissociationConnectivitySet) {
		return;
	}

	/**
	 * The destructor
	 */
	virtual ~Reactant() {
	}

	/**
	 * Note that we result from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param a Number that can be used by daughter classes.
	 * @param b Number that can be used by daughter classes.
	 */
	virtual void resultFrom(ProductionReaction& reaction,
			int a[4] = defaultInit, int b[4] = defaultInit) override {
		return;
	}

	/**
	 * Note that we result from the given reaction involving a super cluster.
	 * Assumes the reaction is already in the network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param prInfos Production reaction parameters used by derived classes.
	 */
	virtual void resultFrom(ProductionReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos)
					override {
		// Must be defined because we use stock Reactants with dummy
		// Reactions, so we need to be able to create Reactant objects.
		return;
	}

	/**
	 * Note that we result from the given reaction involving a super cluster.
	 * Assumes the reaction is already in the network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param product The cluster created by the reaction.
	 */
	virtual void resultFrom(ProductionReaction& reaction, IReactant& product)
			override {
		// Must be defined because we use stock Reactants with dummy
		// Reactions, so we need to be able to create Reactant objects.
		return;
	}

	/**
	 * Note that we result from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param coef Number that can be used by daughter classes.
	 */
	virtual void resultFrom(ProductionReaction& reaction, double *coef)
			override {
		return;
	}

	/**
	 * Note that we combine with another cluster in a production reaction.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param a Number that can be used by daughter classes.
	 */
	virtual void participateIn(ProductionReaction& reaction, int a[4] =
			defaultInit) override {
		return;
	}

	/**
	 * Note that we combine with another cluster in a production reaction
	 * involving a super cluster.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param prInfos Production reaction parameters.
	 */
	virtual void participateIn(ProductionReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos)
					override {
		// Must be defined because we use stock Reactants with dummy
		// Reactions, so we need to be able to create Reactant objects.
		return;
	}

	/**
	 * Note that we combine with another cluster in a production reaction
	 * involving a super cluster.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param product The cluster created by the reaction.
	 */
	virtual void participateIn(ProductionReaction& reaction, IReactant& product)
			override {
		// Must be defined because we use stock Reactants with dummy
		// Reactions, so we need to be able to create Reactant objects.
		return;
	}

	/**
	 * Note that we combine with another cluster in a production reaction.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param coef Number that can be used by daughter classes.
	 */
	virtual void participateIn(ProductionReaction& reaction, double *coef)
			override {
		return;
	}

	/**
	 * Note that we combine with another cluster in a dissociation reaction.
	 * Assumes the reaction is already inour network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param a Number that can be used by daughter classes.
	 * @param b Number that can be used by daughter classes.
	 */
	virtual void participateIn(DissociationReaction& reaction, int a[4] =
			defaultInit, int b[4] = defaultInit) override {
		return;
	}

	/**
	 * Note that we combine with another cluster in a dissociation reaction
	 * involving a super cluster.
	 * Assumes the reaction is already inour network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param prInfos Production reaction parameters.
	 */
	virtual void participateIn(DissociationReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos)
					override {
		// Must be defined because we use stock Reactants with dummy
		// Reactions, so we need to be able to create Reactant objects.
		return;
	}

	/**
	 * Note that we combine with another cluster in a dissociation reaction
	 * involving a super cluster.
	 * Assumes the reaction is already inour network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param disso The dissociating cluster.
	 */
	virtual void participateIn(DissociationReaction& reaction, IReactant& disso)
			override {
		// Must be defined because we use stock Reactants with dummy
		// Reactions, so we need to be able to create Reactant objects.
		return;
	}

	/**
	 * Note that we combine with another cluster in a dissociation reaction.
	 * Assumes the reaction is already inour network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param coef Number that can be used by daughter classes.
	 */
	virtual void participateIn(DissociationReaction& reaction, double *coef)
			override {
		return;
	}

	/**
	 * Note that we emit from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param a Number that can be used by daughter classes.
	 */
	virtual void emitFrom(DissociationReaction& reaction,
			int a[4] = defaultInit) override {
		return;
	}

	/**
	 * Note that we emit from the given reaction involving a super cluster.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param prInfos Production reaction parameters.
	 */
	virtual void emitFrom(DissociationReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos)
					override {
		// Must be defined because we use stock Reactants with dummy
		// Reactions, so we need to be able to create Reactant objects.
		return;
	}

	/**
	 * Note that we emit from the given reaction involving a super cluster.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param disso The dissociating cluster.
	 */
	virtual void emitFrom(DissociationReaction& reaction, IReactant& disso)
			override {
		// Must be defined because we use stock Reactants with dummy
		// Reactions, so we need to be able to create Reactant objects.
		return;
	}

	/**
	 * Note that we emit from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param coef Number that can be used by daughter classes.
	 */
	virtual void emitFrom(DissociationReaction& reaction, double *coef)
			override {
		return;
	}

	/**
	 * Add the reactions to the network lists.
	 */
	virtual void optimizeReactions() override {
		return;
	}

	/**
	 * This operation returns the current concentration.
	 *
	 * @param distA The first distance for super clusters
	 * @param distB The second distance for super clusters
	 * @param distC The third distance for super clusters
	 * @param distD The fourth distance for super clusters
	 * @return The concentration of this reactant
	 */
	virtual double getConcentration(double distA = 0.0, double distB = 0.0,
			double distC = 0.0, double distD = 0.0) const override {
		return concentration;
	}

	/**
	 * This operation sets the concentration of the reactant to the
	 * specified amount.
	 *
	 * @param conc The new concentation
	 */
	void setConcentration(double conc) override {
		concentration = conc;
	}

	/**
	 * This operation returns the total flux of this reactant in the
	 * current network.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The total change in flux for this reactant due to all
	 * reactions
	 */
	virtual double getTotalFlux(int i) override {
		return 0.0;
	}

	/**
	 * Update reactant using other reactants in its network.
	 */
	virtual void updateFromNetwork() override {
		// Nothing to do - derived classes do any meaningful work.
		// Required to be defined because we create explicit Reactant objects,
		// e.g. as dummy objects.
	}

	/**
	 * This operation signifies that the reactant with reactant Id should be
	 * listed as connected with this reactant through forward reactions.
	 *
	 * @param id The integer id of the reactant that is connected
	 * to this reactant
	 */
	void setReactionConnectivity(int id) override {
		reactionConnectivitySet.insert(id);
	}

	/**
	 * This operation signifies that the reactant with reactant Id should be
	 * listed as connected with this reactant through forward reactions.
	 *
	 * @param id The integer id of the reactant that is connected
	 * to this reactant
	 */
	void setDissociationConnectivity(int id) override {
		dissociationConnectivitySet.insert(id);
	}

	/**
	 * This operation reset the connectivity sets based on the information
	 * in the effective production and dissociation vectors.
	 */
	virtual void resetConnectivities() override {
		return;
	}

	/**
	 * Add grid points to the vector of diffusion coefficients or remove
	 * them if the value is negative.
	 *
	 * @param i The number of grid point to add or remove
	 */
	virtual void addGridPoints(int i) override;

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
	virtual std::vector<int> getConnectivity() const override;

	/**
	 * This operation returns the list of partial derivatives of this reactant
	 * with respect to all other reactants in the network. The combined lists
	 * of partial derivatives from all of the reactants in the network can be
	 * used to form, for example, a Jacobian.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return the partial derivatives for this reactant where index zero
	 * corresponds to the first reactant in the list returned by the
	 * ReactionNetwork::getAll() operation.
	 */
	virtual std::vector<double> getPartialDerivatives(int i) const override {
		return std::vector<double>(network.getDOF(), 0.0);
	}

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
	 * @param i The location on the grid in the depth direction
	 */
	virtual void getPartialDerivatives(std::vector<double> & partials,
			int i) const override {
		// nothing to do.
	}

	/**
	 * This operation returns the name of the reactant.
	 *
	 * @return The name
	 */
	const std::string getName() const override {
		return name;
	}

	/**
	 * This operation returns the reactant's type. It is up to subclasses to
	 * define exactly what the allowed types may be.
	 *
	 * @return The type of this reactant.
	 */
	ReactantType getType() const override {
		return type;
	}

	/**
	 * This operation returns the composition of this reactant. This map is empty
	 * when returned by the base class.
	 *
	 * @return The composition returned as a map with keys naming distinct
	 * elements and values indicating the amount of the element present.
	 */
	virtual const IReactant::Composition & getComposition() const override {
		return composition;
	}

	/**
	 * This operation sets the id of the reactant, The id is zero by default
	 * and clients, most likely the ReactionNetwork, are expected to set the
	 * id as needed.
	 *
	 * @param nId The new id for this reactant
	 */
	void setId(int nId) override {
		id = nId;
	}

	/**
	 * This operation returns the id for this reactant.
	 *
	 * @return The id
	 */
	int getId() const override {
		return id;
	}

	/**
	 * This operation sets the id of the first moment of the reactant.
	 *
	 * @param nId The new id for this moment
	 * @param axis The direction
	 */
	void setMomentId(int nId, int axis = 0) override {
		momId[axis] = nId;
	}

	/**
	 * This operation returns the id for this reactant first moment.
	 *
	 * @param axis The direction
	 * @return The id
	 */
	int getMomentId(int axis = 0) const override {
		return momId[axis];
	}

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
	 * @param i The location on the grid
	 */
	void setTemperature(double temp, int i) override;

	/**
	 * This operation returns the temperature at which the reactant currently exists.
	 *
	 * @param i The location on the grid
	 * @return The temperature.
	 */
	double getTemperature(int i) const override {
		return temperature[i];
	}

	/**
	 * This operation returns the total size of the reactant.
	 *
	 * @return The total size of this reactant including the contributions
	 * from all species types
	 */
	IReactant::SizeType getSize() const override {
		return size;
	}

	/**
	 * This operation retrieves the formation energy for this reactant.
	 *
	 * @return The value of the formation energy
	 */
	double getFormationEnergy() const override {
		return formationEnergy;
	}

	/**
	 * This operation sets the formation energy for this reactant.
	 *
	 * @param energy The formation energy
	 */
	void setFormationEnergy(double energy) override {
		formationEnergy = energy;
	}

	/**
	 * This operation retrieves the diffusion factor, D_0, that is used to
	 * calculate the diffusion coefficient for this reactant.
	 *
	 * @return The diffusion factor of this reactant
	 */
	double getDiffusionFactor() const override {
		return diffusionFactor;
	}

	/**
	 * This operation sets the diffusion factor, D_0, that is used to calculate
	 * the diffusion coefficient for this reactant.
	 *
	 * @param factor The diffusion factor
	 */
	virtual void setDiffusionFactor(const double factor) override;

	/**
	 * This operation returns the diffusion coefficient for this reactant and is
	 * calculated from the diffusion factor.
	 *
	 * @param i The position on the grid
	 * @return The diffusion coefficient
	 */
	double getDiffusionCoefficient(int i) const override {
		return diffusionCoefficient[i];
	}

	/**
	 * This operation sets the migration energy for this reactant.
	 *
	 * @param energy The migration energy
	 */
	virtual void setMigrationEnergy(const double energy) override;

	/**
	 * This operation retrieves the migration energy for this reactant.
	 *
	 * @return the migration energy
	 */
	double getMigrationEnergy() const override {
		return migrationEnergy;
	}

	/**
	 * This operation returns the reaction radius for the
	 * particular reactant.
	 *
	 * @return The reaction radius
	 */
	double getReactionRadius() const override {
		return reactionRadius;
	}

	/**
	 * This operation returns the sum of combination rate and emission rate
	 * (where this reactant is on the left side of the reaction) for this
	 * particular reactant.
	 * This is used to computed the desorption rate in the
	 * modified trap-mutation handler.
	 *
	 * @param i The position on the grid
	 * @return The rate
	 */
	virtual double getLeftSideRate(int i) const override {
		return 0.0;
	}

	/**
	 * This operation returns the vector of production reactions in which
	 * this cluster is involved, containing the id of the reactants, the rate, and
	 * the coefs[0][0]
	 *
	 * @return The vector of productions
	 */
	virtual std::vector<std::vector<double> > getProdVector() const override {
		return std::vector<std::vector<double> >();
	}

	/**
	 * This operation returns the vector of combination reactions in which
	 * this cluster is involved, containing the id of the other reactants, the rate, and
	 * the coefs[0]
	 *
	 * @return The vector of combinations
	 */
	virtual std::vector<std::vector<double> > getCombVector() const override {
		return std::vector<std::vector<double> >();
	}

	/**
	 * This operation returns the vector of dissociation reactions in which
	 * this cluster is involved, containing the id of the emitting reactants, the rate, and
	 * the coefs[0][0]
	 *
	 * @return The vector of dissociations
	 */
	virtual std::vector<std::vector<double> > getDissoVector() const override {
		return std::vector<std::vector<double> >();
	}

	/**
	 * This operation returns the vector of emission reactions in which
	 * this cluster is involved, containing the rate, and
	 * the coefs[0][0]
	 *
	 * @return The vector of productions
	 */
	virtual std::vector<std::vector<double> > getEmitVector() const override {
		return std::vector<std::vector<double> >();
	}

	/**
	 * This operation returns true if the cluster is a mixed-species or compound
	 * cluster and false if it is a single species cluster.
	 */
	virtual bool isMixed() const override {
		return false;
	}

	/**
	 * Tell reactant to output a representation of its reaction coefficients
	 * to the given output stream.
	 *
	 * @param os Output stream on which to output coefficients.
	 */
	// We must define this because the code may use a stock Reactant
	// when using dummy reactions, and thus we have to define all
	// pure virtual functions from our base class(es).
	virtual void outputCoefficientsTo(std::ostream& os) const override {
		// Nothing to do.
	}
};

} // end namespace xolotlCore

#endif
