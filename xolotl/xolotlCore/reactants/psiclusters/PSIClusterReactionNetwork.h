#ifndef PSI_CLUSTER_REACTION_NETWORK_H
#define PSI_CLUSTER_REACTION_NETWORK_H

// Includes
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <algorithm>
#include "ReactionNetwork.h"
#include "PSISuperCluster.h"
#include "ReactantType.h"

namespace xolotlCore {

/**
 *  This class manages the set of reactants and compound reactants (
 *  combinations of normal reactants) for PSI clusters. It also manages a
 *  set of properties that describes the total collection.
 *
 *  This class is a very heavyweight class that should not be abused.
 *
 *  Reactants that are added to this network must be added as with shared_ptrs.
 *  Furthermore, reactants that are added to this network have their ids set to
 *  a network specific id. Reactants should not be shared between separate
 *  instances of a PSIClusterReactionNetwork.
 */
class PSIClusterReactionNetwork: public ReactionNetwork {

private:
	/**
	 * Concise name for map supporting quick lookup of supercluster containing
	 * specifc number of He and V.
	 *
	 * We could use a map, but because we expect it to be dense (i.e.,
	 * most pairs of He and V counts will have a valid super cluster),
	 * a 2D matrix indexed by nHe and nV gives better performance for
	 * lookups without costing too much more (if any) in terms of memory.
	 * We use a vector of vectors instead of array of arrays because
	 * we don't know the dimensions at compile time.  We use a vector of
	 * vectors instead of a single vector with more complex indexing
	 * to simplify indexing.  It may be worthwhile to compare performance
	 * with a single vector implementation to see if the single
	 * memory allocation ends up giving better performance.
	 */
	using HeVToSuperClusterMap = std::vector<std::vector<ReactantMap::const_iterator> >;

	/**
	 * Map supporting quick identification of super cluster containing
	 * given number of He and V.
	 */
	HeVToSuperClusterMap superClusterLookupMap;

	/**
	 * Calculate the dissociation constant of the first cluster with respect to
	 * the single-species cluster of the same type based on the current clusters
	 * atomic volume, reaction rate constant, and binding energies.
	 *
	 * @param reaction The reaction
	 * @return The dissociation constant
	 */
	double calculateDissociationConstant(
			const DissociationReaction& reaction) const override;

	/**
	 * Calculate the binding energy for the dissociation cluster to emit the single
	 * and second cluster.
	 *
	 * @param reaction The reaction
	 * @return The binding energy corresponding to this dissociation
	 */
	double computeBindingEnergy(const DissociationReaction& reaction) const
			override;

	/**
	 * Find the super cluster that contains the original cluster with nHe
	 * helium atoms and nV vacancies.
	 *
	 * @param nHe The number of helium atoms
	 * @param nD The number of deuterium atoms
	 * @param nT The number of tritium atoms
	 * @param nV The number of vacancies
	 * @return The super cluster representing the cluster with nHe helium
	 * and nV vacancies, or nullptr if no such cluster exists.
	 */
	IReactant * getSuperFromComp(IReactant::SizeType nHe,
			IReactant::SizeType nD, IReactant::SizeType nT,
			IReactant::SizeType nV);

	ProductionReaction& defineReactionBase(IReactant& r1, IReactant& r2,
			int a[4] = defaultInit, bool secondProduct = false)
					__attribute__((always_inline)) {

		// Add a production reaction to our network.
		std::unique_ptr<ProductionReaction> reaction(
				new ProductionReaction(r1, r2));
		auto& prref = add(std::move(reaction));

		// Tell the reactants that they are involved in this reaction
		// if this is not the second product
		if (secondProduct)
			return prref;

		r1.participateIn(prref, a);
		r2.participateIn(prref, a);

		return prref;
	}

	void defineAnnihilationReaction(IReactant& r1, IReactant& r2,
			IReactant& product) {

		// Define the basic reaction.
		auto& reaction = defineReactionBase(r1, r2);

		// Tell the product it results from the reaction.
		product.resultFrom(reaction);
	}

	void defineCompleteAnnihilationReaction(IReactant& r1, IReactant& r2) {

		// Define the basic reaction
		defineReactionBase(r1, r2);

		// Nothing else to do since there is no product.
	}

	void defineProductionReaction(IReactant& r1, IReactant& super,
			IReactant& product, int a[4] = defaultInit, int b[4] = defaultInit,
			bool secondProduct = false) {
		// Define the basic production reaction.
		auto& reaction = defineReactionBase(r1, super, b, secondProduct);

		// Tell product it is a product of this reaction.
		product.resultFrom(reaction, a, b);

		if (secondProduct)
			return;

		// Check if reverse reaction is allowed.
		checkForDissociation(product, reaction, a, b);
	}

	/**
	 * Define a batch of production reactions for the given
	 * pair of reactants.
	 *
	 * @param r1 A reactant involved in a production reaction.
	 * @param r2 The super reactant involved in a production reaction.
	 * @param pris Information about reactants are involved with each reaction.
	 * @param secondProduct If we are setting the reaction for the second product.
	 */
	void defineProductionReactions(IReactant& r1, IReactant& super,
			const std::vector<PendingProductionReactionInfo>& pris,
			bool secondProduct = false);

	// TODO should we default a, b, c, d to 0?
	void defineDissociationReaction(ProductionReaction& forwardReaction,
			IReactant& emitting, int a[4] = defaultInit,
			int b[4] = defaultInit) {

		std::unique_ptr<DissociationReaction> dissociationReaction(
				new DissociationReaction(emitting, forwardReaction.first,
						forwardReaction.second, &forwardReaction));
		auto& drref = add(std::move(dissociationReaction));

		// Tell the reactants that they are in this reaction
		forwardReaction.first.participateIn(drref, a, b);
		forwardReaction.second.participateIn(drref, a, b);
		emitting.emitFrom(drref, a);
	}

	/**
	 * Define a batch of dissociation reactions for the given
	 * forward reaction.
	 *
	 * @param forwardReaction The forward reaction in question.
	 * @param prodMap Map of reaction parameters, keyed by the product
	 * of the reaction.
	 */
	// TODO possible to use a ref for the key?
	using ProductToProductionMap =
	std::unordered_map<IReactant*, std::vector<PendingProductionReactionInfo> >;

	void defineDissociationReactions(ProductionReaction& forwardReaction,
			const ProductToProductionMap& prodMap);

	/**
	 * Check whether dissociation reaction is allowed for
	 * given production reaction.
	 *
	 * @param emittingReactant The reactant that would emit the pair
	 * @param reaction The reaction to test.
	 * @return true iff dissociation for the given reaction is allowed.
	 */
	bool canDissociate(IReactant& emittingReactant,
			ProductionReaction& reaction) const;

	/**
	 * Add the dissociation connectivity for the reverse reaction if it is allowed.
	 *
	 * @param emittingReactant The reactant that would emit the pair
	 * @param reaction The reaction we want to reverse
	 * @param a The helium number for the emitting superCluster
	 * @param b The vacancy number for the emitting superCluster
	 *
	 */
	void checkForDissociation(IReactant& emittingReactant,
			ProductionReaction& reaction, int a[4] = defaultInit, int b[4] =
					defaultInit);

public:

	/**
	 * Default constructor, deleted to force construction using parameters.
	 */
	PSIClusterReactionNetwork() = delete;

	/**
	 * The Constructor
	 *
	 * @param registry The performance handler registry
	 */
	PSIClusterReactionNetwork(
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	PSIClusterReactionNetwork(const PSIClusterReactionNetwork& other) = delete;

	/**
	 * Computes the full reaction connectivity matrix for this network.
	 */
	void createReactionConnectivity();

	/**
	 * This operation sets the temperature at which the reactants currently
	 * exists. It calls setTemperature() on each reactant.
	 *
	 * This is the simplest way to set the temperature for all reactants is to
	 * call the ReactionNetwork::setTemperature() operation.
	 *
	 * @param temp The new temperature
	 */
	void setTemperature(double temp) override;

	/**
	 * This operation reinitializes the network.
	 *
	 * It computes the cluster Ids and network size from the allReactants vector.
	 */
	void reinitializeNetwork() override;

	/**
	 * This method redefines the connectivities for each cluster in the
	 * allReactans vector.
	 */
	void reinitializeConnectivities() override;

	/**
	 * This operation updates the concentrations for all reactants in the
	 * network from an array.
	 *
	 * @param concentrations The array of doubles that will be for the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. Properly aligning the array in memory so that this operation
	 * does not overrun is up to the caller.
	 */
	void updateConcentrationsFromArray(double * concentrations) override;

	/**
	 * This operation returns the size or number of reactants and momentums in the network.
	 *
	 * @return The number of degrees of freedom
	 */
	virtual int getDOF() const override {
		return size() + 3 * getAll(ReactantType::PSISuper).size() + 1;
	}

	/**
	 * This operation returns the list (vector) of each reactant in the network.
	 *
	 * @return The list of compositions
	 */
	virtual std::vector<std::vector<int> > getCompositionList() const override;

	/**
	 * Get the diagonal fill for the Jacobian, corresponding to the reactions.
	 *
	 * @param diagFill The pointer to the vector where the connectivity information is kept
	 */
	void getDiagonalFill(int *diagFill) override;

	/**
	 * Get the total concentration of atoms contained in the network.
	 *
	 * Here the atoms that are considered are:
	 * 0 helium
	 * 1 deuterium
	 * 2 tritium
	 *
	 * @param i Index to switch between the different types of atoms
	 * @return The total concentration
	 */
	double getTotalAtomConcentration(int i = 0) override;

	/**
	 * Get the total concentration of atoms contained in bubbles in the network.
	 *
	 * Here the atoms that are considered are:
	 * 0 helium
	 * 1 deuterium
	 * 2 tritium
	 *
	 * @param i Index to switch between the different types of atoms
	 * @return The total concentration
	 */
	double getTotalTrappedAtomConcentration(int i = 0) override;

	/**
	 * Get the total concentration of vacancies contained in the network.
	 *
	 * @return The total concentration
	 */
	double getTotalVConcentration() override;

	/**
	 * Get the total concentration of tungsten interstitials in the network.
	 *
	 * @return The total concentration
	 */
	double getTotalIConcentration() override;

	/**
	 * Calculate all the rate constants for the reactions and dissociations of the network.
	 * Need to be called only when the temperature changes.
	 */
	void computeRateConstants() override;

	/**
	 * Compute the fluxes generated by all the reactions
	 * for all the clusters and their momentums.
	 *
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the fluxes are computed used to find the next solution
	 */
	void computeAllFluxes(double *updatedConcOffset) override;

	/**
	 * Compute the partial derivatives generated by all the reactions
	 * for all the clusters and their momentum.
	 *
	 * @param vals The pointer to the array that will contain the values of
	 * partials for the reactions
	 * @param indices The pointer to the array that will contain the indices
	 * of the clusters
	 * @param size The pointer to the array that will contain the number of reactions for
	 * this cluster
	 */
	virtual void computeAllPartials(double *vals, int *indices, int *size)
			override;
};

} // namespace xolotlCore

#endif
