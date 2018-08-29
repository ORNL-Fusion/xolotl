#ifndef REACTION_NETWORK_H
#define REACTION_NETWORK_H

// Includes
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <Constants.h>
#include "IReactionNetwork.h"
#include "Reactant.h"

namespace xolotlPerf {
class IHandlerRegistry;
class IEventCounter;
}

// Using std::unordered_map often gives better performance than std::map,
// but requires us to define our own hash function for more complex types.
using ReactantSizePair = std::pair<xolotlCore::IReactant::SizeType, xolotlCore::IReactant::SizeType>;

namespace std {

template<>
struct hash<ReactantSizePair> {
	size_t operator()(const ReactantSizePair& pr) const {
		return (pr.first << 16) + (pr.second);
	}
};

} // namespace std

namespace xolotlCore {

/**
 *  This class manages the set of reactants and compound reactants
 *  (combinations of normal reactants). It also manages a set of properties
 *  that describe both.
 */
class ReactionNetwork: public IReactionNetwork {

public:
	using PartialsIdxMap = std::unordered_map<int, int>;

private:
	/**
	 * Types of reactants that we support.
	 */
	std::set<ReactantType> knownReactantTypes;

protected:

	/**
	 * A functor useful for identifying a set of reactants by their
	 * composition from a container, e.g., when removing a collection
	 * of doomed reactants from a vector of reactants.
	 */
	class ReactantMatcher {
	private:
		/**
		 * The canonical composition representations of the reactants
		 * we are to find.
		 */
		std::set<IReactant::Composition> comps;

	public:
		/**
		 * Build a ReactantMatcher.
		 * @param reactants The collection of reactants we are to recognize.
		 */
		ReactantMatcher(const IReactionNetwork::ReactantVector& reactants) {
			for (IReactant const& reactant : reactants) {
				comps.insert(reactant.getComposition());
			}
		}

		/**
		 * Determine if the given reactant is in our set.
		 * @param testReactant The reactant to check.
		 * @return true iff the reactant's composition's canonical
		 * representation is in our set.
		 */
		bool operator()(const IReactant* testReactant) const {
			auto iter = comps.find(testReactant->getComposition());
			return (iter != comps.end());
		}

		/**
		 * Determine if the given reactant is in our set.
		 * @param testReactant The reactant to check.
		 * @return true iff the reactant's composition's canonical
		 * representation is in our set.
		 */
		bool operator()(const IReactant& testReactant) const {
			auto iter = comps.find(testReactant.getComposition());
			return (iter != comps.end());
		}

		/**
		 * Determine if the given reactant is in our set.
		 * @param testReactant The reactant to check.
		 * @return true iff the reactant's composition's canonical 
		 * representation is in our set of doomed reactants.
		 */
		bool operator()(const std::unique_ptr<IReactant>& testReactant) const {
			return this->operator()(testReactant.get());
		}

		/**
		 * Determine if the given reactant is in our set.
		 * @param testMapItem ReactantMap item for reactant to check.
		 * @return true iff the reactant composition's representation
		 * is in our set of doomed reactants.
		 */
		bool operator()(const ReactantMap::value_type& testMapItem) const {
			auto const& testReactant = testMapItem.second;
			return this->operator()(*testReactant);
		}
	};

	/**
	 * The performance handler registry that will be used with
	 * this class.
	 */
	std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry;

	/**
	 * All known ProductionReactions in the network, keyed by a
	 * representation of the reaction.
	 */
	using ProductionReactionMap = std::unordered_map<ProductionReaction::KeyType, std::unique_ptr<ProductionReaction> >;
	ProductionReactionMap productionReactionMap;

	/**
	 * All known dissociation reactions in the network, keyed by a
	 * representation of the reaction.
	 */
	using DissociationReactionMap = std::unordered_map<DissociationReaction::KeyType,
	std::unique_ptr<DissociationReaction> >;
	DissociationReactionMap dissociationReactionMap;

	/**
	 * A map for storing the dfill configuration and accelerating the formation of
	 * the Jacobian. Its keys are reactant/cluster ids and its values are integer
	 * vectors of the column ids that are marked as connected for that cluster in
	 * the dfill array.
	 */
	std::unordered_map<int, std::vector<int> > dFillMap;

	/**
	 * Map of reactant/cluster ids to another map of ids in the DOF space
	 * to its location within the sparse partials arrays.
	 */
	std::unordered_map<int, PartialsIdxMap> dFillInvMap;

	/**
	 * The current temperature at which the network's clusters exist.
	 */
	double temperature;

	/**
	 * The biggest rate for this cluster
	 */
	double biggestRate;

	/**
	 * Are dissociations enabled?
	 */
	bool dissociationsEnabled;

	/**
	 * Maximum cluster sizes currently in the network for each
	 * of the reactant types we support.
	 */
	std::map<ReactantType, IReactant::SizeType> maxClusterSizeMap;

	/**
	 * This vector contains the information on the group bounds in both directions.
	 */
	std::vector<IReactant::SizeType> boundVector;

	/**
	 * All reactants known to the network.
	 */
	IReactant::RefVector allReactants;

	/**
	 * This map stores all of the clusters in the network by type.
	 */
	std::unordered_map<ReactantType, ReactantMap> clusterTypeMap;

	/**
	 * Calculate the reaction constant dependent on the
	 * reaction radii and the diffusion coefficients for the
	 * ith and jth clusters, which itself depends on the current
	 * temperature.
	 *
	 * @param reaction The reaction
	 * @param i The location on the grid
	 * @return The rate
	 */
	virtual double calculateReactionRateConstant(const ProductionReaction& reaction,
			int i) const;

	/**
	 * Calculate the dissociation constant of the first cluster with respect to
	 * the single-species cluster of the same type based on the current clusters
	 * atomic volume, reaction rate constant, and binding energies.
	 *
	 * Need to be overwritten by daughter classes because of the atomic volume.
	 *
	 * @param reaction The reaction
	 * @param i The location on the grid in the depth direction
	 * @return The dissociation constant
	 */
	virtual double calculateDissociationConstant(
			const DissociationReaction& reaction, int i) = 0;

	/**
	 * Calculate the binding energy for the dissociation cluster to emit the single
	 * and second cluster.
	 *
	 * @param reaction The reaction
	 * @return The binding energy corresponding to this dissociation
	 */
	virtual double computeBindingEnergy(
			const DissociationReaction& reaction) const = 0;

	/**
	 * Find index of interval in boundVector that contains a value.
	 * Assumes that:
	 *   boundVector is sorted
	 *   boundVector indicates non-overlapping intervals
	 *   boundVector last element indicates the upper bound of the last
	 interval (one past, following usual C++ standard library convention.)
	 *   0 is not included in the first boundVector interval.
	 *
	 * @param val  The value of interest.
	 * @return The lower bound of the interval in boundVector that contains
	 *       'val.'  If there is no such interval, returns 0.
	 */
	std::size_t findBoundsIntervalBaseIdx(IReactant::SizeType val) const {

		// Find the first item that is *greater than* the given count.
		auto iter = std::upper_bound(boundVector.begin(), boundVector.end(),
				val);

		// There are three cases for iter:
		// * std::upper_bound returned begin() => val is smaller than any interval.
		// * std::upper_bound returned end() => val is larger than any interval.
		// * otherwise => iter points to *next* item after interval we want.
		return ((iter != boundVector.begin()) and (iter != boundVector.end())) ?
				(iter - boundVector.begin()) :
				std::numeric_limits<std::size_t>::max();
	}

public:

	/**
	 * Default constructor, deleted because we need params to construct.
	 */
	ReactionNetwork() = delete;

	/**
	 * The constructor that takes the performance handler registry.
	 * It initializes the properties and reactants vector.
	 *
	 * @param _knownReactantTypes Reactant types that we support.
	 * @param _registry The performance handler registry
	 */
	ReactionNetwork(const std::set<ReactantType>& _knownReactantTypes,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> _registry);

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	ReactionNetwork(const ReactionNetwork& other) = delete;

	/**
	 * The destructor.
	 */
	virtual ~ReactionNetwork() {
	}

	/**
	 * This operation sets the temperature at which the reactants currently
	 * exists. It calls setTemperature() on each reactant.
	 *
	 * This is the simplest way to set the temperature for all reactants.
	 *
	 * @param temp The new temperature
	 * @param i The location on the grid in the depth direction
	 */
	virtual void setTemperature(double temp, int i = 0) override;

	/**
	 * This operation returns the temperature at which the cluster currently exists.
	 *
	 * @return The temperature
	 */
	virtual double getTemperature() const override {
		return temperature;
	}

	/**
	 * Retrieve the single-species reactant with the given type and size if it
	 * exists in the network or null if not.
	 * Convenience function for get() that takes a
	 * reactant type and composition.
	 *
	 * @param species The reactant's single species.
	 * @param size The size of the reactant.
	 * @return A pointer to the reactant, or nullptr if it does not 
	 * exist in the network.
	 */
	IReactant * get(Species species, IReactant::SizeType size) const override {
		IReactant::Composition comp;
		comp[toCompIdx(species)] = size;
		return get(toReactantType(species), comp);
	}

	/**
	 * Retrieve the reactant with the given type and composition if
	 * exists in the network.
	 *
	 * @param type The type of the reactant
	 * @param comp The composition of the reactant
	 * @return A pointer to the reactant of type 'type' and with composition
	 * 'comp.' nullptr if no such reactant exists.
	 */
	IReactant * get(ReactantType type, const IReactant::Composition& comp) const
			override;

	/**
	 * This operation returns all reactants in the network without regard for
	 * their composition or whether they are compound reactants. The list may
	 * or may not be ordered and the decision is left to implementers.
	 *
	 * @return The list of all of the reactants in the network
	 */
	virtual const IReactant::RefVector& getAll() const override {

		return allReactants;
	}

	/**
	 * This operation returns all reactants in the network with the given type.
	 * The list may or may not be ordered and the decision is left to
	 * implementers.
	 *
	 * @param type The reactant or compound reactant type
	 * @return The list of all of the reactants in the network or null if the
	 * type is invalid
	 */
	virtual const ReactantMap& getAll(ReactantType type) const override {

		return clusterTypeMap.at(type);
	}

	/**
	 * Give the reactant to the network.
	 *
	 * This operation sets the id of the reactant to one that is specific
	 * to this network. Do not share reactants across networks! This id is
	 * guaranteed to be between 1 and n, including both, for n reactants in
	 * the network.
	 *
	 * @param reactant The reactant that should be added to the network.
	 */
	void add(std::unique_ptr<IReactant> reactant) override;

	/**
	 * This operation reinitializes the network.
	 *
	 * It computes the cluster Ids and network size from the allReactants vector.
	 */
	virtual void reinitializeNetwork() override {
		return;
	}

	/**
	 * This operation returns the size or number of reactants in the network.
	 *
	 * @return The number of reactants in the network
	 */
	int size() const override {
		return allReactants.size();
	}

	/**
	 * This operation returns the number of super reactants in the network.
	 * Needs to be overridden by the daughter classes.
	 *
	 * @return The number of super reactants in the network
	 */
	int getSuperSize() const override {
		return 0;
	}

	/**
	 * This operation returns the size or number of reactants and momentums in the network.
	 *
	 * @return The number of degrees of freedom
	 */
	virtual int getDOF() const override {
		return size();
	}

	/**
	 * This operation returns the list (vector) of each reactant in the network.
	 * Need to be implemented by the daughter classes,
	 *
	 * @return The list of compositions
	 */
	virtual std::vector<std::vector<int> > getCompositionList() const override {
		return std::vector<std::vector<int> >();
	}

	/**
	 * Find the super cluster that contains the original cluster with the given composition.
	 *
	 * @param nHe The number of helium/xenon atoms
	 * @param nD The number of deuterium atoms
	 * @param nT The number of tritium atoms
	 * @param nV The number of vacancies
	 * @return The super cluster representing the cluster with nHe helium
	 * and nV vacancies, or nullptr if no such cluster exists.
	 */
	virtual IReactant * getSuperFromComp(IReactant::SizeType nHe,
			IReactant::SizeType nD, IReactant::SizeType nT,
			IReactant::SizeType nV) const override {
		return nullptr;
	}

	/**
	 * Add a production reaction to the network.
	 *
	 * @param reaction The reaction that should be added to the network
	 * @return The reaction that is now in the network
	 */
	ProductionReaction& add(std::unique_ptr<ProductionReaction> reaction)
			override;

	/**
	 * Add a dissociation reaction to the network.
	 *
	 * @param reaction The reaction that should be added to the network
	 * @return The reaction that is now in the network
	 */
	DissociationReaction& add(std::unique_ptr<DissociationReaction> reaction)
			override;

	/**
	 * This operation fills an array of doubles with the concentrations of all
	 * of the reactants in the network.
	 *
	 * @param concentrations The array that will be filled with the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. If the array is too small to hold the concentrations, SIGSEGV will
	 * be thrown.
	 */
	void fillConcentrationsArray(double * concentrations) override;

	/**
	 * This operation updates the concentrations for all reactants in the
	 * network from an array.
	 *
	 * @param concentrations The array of doubles that will be for the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. Properly aligning the array in memory so that this operation
	 * does not overrun is up to the caller.
	 */
	virtual void updateConcentrationsFromArray(double * concentrations)
			override;

	/**
	 * Get the diagonal fill for the Jacobian, corresponding to the reactions.
	 *
	 * Do nothing here, this method need to be implemented in
	 * subclasses.
	 *
	 * @param sfm Connectivity map.
	 */
	virtual void getDiagonalFill(SparseFillMap& sfm) override {
		return;
	}

	/**
	 * Get the total concentration of atoms contained in the network.
	 *
	 * Returns 0.0 here and needs to be implemented by the daughter classes.
	 *
	 * @param i Index to switch between the different types of atoms
	 * @return The total concentration
	 */
	virtual double getTotalAtomConcentration(int i = 0) override {
		return 0.0;
	}

	/**
	 * Get the total concentration of atoms contained in bubbles in the network.
	 *
	 * Returns 0.0 here and needs to be implemented by the daughter classes.
	 *
	 * @param i Index to switch between the different types of atoms
	 * @return The total concentration
	 */
	virtual double getTotalTrappedAtomConcentration(int i = 0) override {
		return 0.0;
	}

	/**
	 * Get the total concentration of vacancies contained in the network.
	 *
	 * Returns 0.0 here and needs to be implemented by the daughter classes.
	 *
	 * @return The total concentration
	 */
	virtual double getTotalVConcentration() override {
		return 0.0;
	}

	/**
	 * Get the total concentration of material interstitials in the network.
	 *
	 * Returns 0.0 here and needs to be implemented by the daughter classes.
	 *
	 * @return The total concentration
	 */
	virtual double getTotalIConcentration() override {
		return 0.0;
	}

	/**
	 * Calculate all the rate constants for the reactions and dissociations of the network.
	 * Need to be called only when the temperature changes.
	 *
	 * Need to be overwritten by daughter classes.
	 *
	 * @param i The location on the grid in the depth direction
	 */
	virtual void computeRateConstants(int i) override;

	/**
	 * Add grid points to the vector of rates or remove them if the value is negative.
	 *
	 * @param i The number of grid point to add or remove
	 */
	virtual void addGridPoints(int i) override;

	/**
	 * Compute the fluxes generated by all the reactions
	 * for all the clusters and their momentum.
	 *
	 * Do nothing here, this method need to be implemented in
	 * subclasses.
	 *
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the fluxes are computed used to find the next solution
	 * @param i The location on the grid in the depth direction
	 */
	virtual void computeAllFluxes(double *updatedConcOffset, int i = 0)
			override {
		return;
	}

	/**
	 * This operation returns the biggest production rate in the network.
	 *
	 * @return The biggest rate
	 */
	double getBiggestRate() const override {
		return biggestRate;
	}

	/**
	 * Are dissociations enabled?
	 * @return true if reactions are enabled, false otherwise.
	 */
	bool getDissociationsEnabled() const override {
		return dissociationsEnabled;
	}

	/**
	 * This operation returns the phase space list needed to set up the grouping
	 * correctly in PSI.
	 *
	 * @return The phase space list
	 */
	virtual Array<int, 5> getPhaseSpaceList() const override {
		Array<int, 5> list;
		list.Init(0);
		return list;
	}

	/**
	 * Enable dissociations.
	 */
	void enableDissociations() {
		dissociationsEnabled = true;
	}

	/**
	 * Disable dissociations.
	 */
	void disableDissociations() {
		dissociationsEnabled = false;
	}

	/**
	 * Find maximum cluster size currently in the network
	 * for the given reactant type.
	 *
	 * @param rtype Reactant type of interest.
	 * @return Maximum size of cluster of type rtype currently in network.
	 */
	IReactant::SizeType getMaxClusterSize(ReactantType rtype) const override {
		return maxClusterSizeMap.at(rtype);
	}

	/**
	 * Remove the given reactants from the network.
	 *
	 * @param reactants The reactants that should be removed.
	 */
	void removeReactants(const ReactantVector& reactants) override;

	/**
	 * Dump a representation of the network to the given output stream.
	 *
	 * @param os Output stream on which to write network description.
	 */
	void dumpTo(std::ostream& os) const override;

	/**
	 * Obtain reactant types supported by our network.
	 * A reactant type might be returned even though no reactant
	 * of that type currently exists in the network.
	 *
	 * @return Collection of reactant types supported by our network.
	 */
	const std::set<ReactantType>& getKnownReactantTypes() const override {
		return knownReactantTypes;
	}

	/**
	 * Determine the number of partials for each cluster
	 * and their starting locations within the vectors used
	 * when computing partials.
	 *
	 * @param size Number of partials for each cluster.
	 * @param startingIdx Starting index of items owned by each reactant
	 *      within the partials values array and the indices array.
	 * @param i The location on the grid in the depth direction
	 * @return Total number of partials for all clusters.
	 */
	size_t initPartialsSizes(std::vector<int>& size,
			std::vector<size_t>& startingIdx) const override;

	/**
	 * Initialize the indexing used for computing partial derivatives
	 * for all clusters and their momenta.
	 *
	 * @param size Number of partials for each cluster.
	 * @param startingIdx Starting index of items owned by each reactant
	 *      within the partials values array and the indices array.
	 * @param indices The indices of the clusters for the partial derivatives.
	 */
	void initPartialsIndices(const std::vector<int>& size,
			const std::vector<size_t>& startingIdx,
			std::vector<int>& indices) const override;

};

}

#endif
