#ifndef URZ_CLUSTER_REACTION_NETWORK_H
#define URZ_CLUSTER_REACTION_NETWORK_H

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
#include "ReactantType.h"

namespace xolotlCore {

/**
 *  This class manages the set of reactants and compound reactants (
 *  combinations of normal reactants) for UZr clusters. It also manages a
 *  set of properties that describes the total collection.
 *
 *  This class is a very heavyweight class that should not be abused.
 *
 *  Reactants that are added to this network must be added as with shared_ptrs.
 *  Furthermore, reactants that are added to this network have their ids set to
 *  a network specific id. Reactants should not be shared between separate
 *  instances of a UZrClusterReactionNetwork.
 */
class UZrClusterReactionNetwork: public ReactionNetwork {

private:

	//! The fission rate in nm2/s
	double fissionRate;

	/**
	 * Calculate the dissociation constant of the first cluster with respect to
	 * the single-species cluster of the same type based on the current clusters
	 * atomic volume, reaction rate constant, and binding energies.
	 *
	 * @param reaction The reaction
	 * @param i The location on the grid in the depth direction
	 * @return The dissociation constant
	 */
	double calculateDissociationConstant(
			const DissociationReaction& reaction, int i) override;

	/**
	 * Calculate the binding energy for the dissociation cluster to emit the single
	 * and second cluster.
	 *
	 * @param reaction The reaction
	 * @return The binding energy corresponding to this dissociation
	 */
	double computeBindingEnergy(const DissociationReaction& reaction) const
			override;

	void defineProductionReaction(IReactant& r1, IReactant& r2,
			IReactant& product) {
		// Create the reaction
		std::unique_ptr<ProductionReaction> reaction(
				new ProductionReaction(r1,
						r2));
		auto& prref = add(std::move(reaction));
		// Tell the reactants that they are in this reaction
		r1.participateIn(prref, product);
		r2.participateIn(prref, product);
		product.resultFrom(prref, product);

		// Check for dissociation
		if (r1.getSize() == 1 || r2.getSize() == 1) {
			defineDissociationReaction(prref, product);
		}
	}

	void defineDissociationReaction(ProductionReaction& forwardReaction,
			IReactant& emitting) {
		// Create the reaction
		std::unique_ptr<DissociationReaction> dissociationReaction(
				new DissociationReaction(emitting, forwardReaction.first,
						forwardReaction.second));
		// Set the reverse reaction
		dissociationReaction->reverseReaction = &forwardReaction;
		auto& drref = add(std::move(dissociationReaction));
		// Tell the reactants that their are in this reaction
		forwardReaction.first.participateIn(drref, emitting);
		forwardReaction.second.participateIn(drref, emitting);
		emitting.emitFrom(drref, emitting);
	}

public:

	/**
	 * Default constructor, deleted to force construction using parameters.
	 */
	UZrClusterReactionNetwork() = delete;

	/**
	 * The Constructor
	 *
	 * @param registry The performance handler registry
	 */
	UZrClusterReactionNetwork(
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	UZrClusterReactionNetwork(const UZrClusterReactionNetwork& other) = delete;

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
	 * @param i The location on the grid in the depth direction
	 */
	virtual void setTemperature(double temp, int i) override;

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
		return size() + 1;
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
	 * @param sfm Connectivity map.
	 */
	void getDiagonalFill(SparseFillMap& sfm) override;

	/**
	 * Get the total concentration of atoms contained in the network.
	 *
	 * Here the atoms that are considered are helium atoms.
	 *
	 * @param i Index to switch between the different types of atoms
	 * @return The total concentration
	 */
	double getTotalAtomConcentration(int i = 0) override;

	/**
	 * Get the total concentration of atoms contained in bubbles in the network.
	 *
	 * Here the atoms that are considered are helium atoms.
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
	 * This operation sets the fission rate, needed to compute the diffusion coefficient
	 * in NE.
	 *
	 * @param rate The fission rate
	 */
	void setFissionRate(double rate) override {
		fissionRate = rate;
		return;
	}

	/**
	 * This operation returns the fission rate, needed to compute the diffusion coefficient
	 * in NE.
	 *
	 * @return The fission rate
	 */
	double getFissionRate() const override {
		return fissionRate;
	}

	/**
	 * Compute the fluxes generated by all the reactions
	 * for all the clusters and their momentums.
	 *
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the fluxes are computed used to find the next solution
	 * @param i The location on the grid in the depth direction
	 */
	void computeAllFluxes(double *updatedConcOffset, int i) override;

	/**
	 * Compute the partial derivatives generated by all the reactions
	 * for all the clusters and their momentum.
	 *
	 * @param startingIdx Starting index of items owned by each reactant
	 *      within the partials values array and the indices array.
	 * @param indices The indices of the clusters for the partial derivatives.
	 * @param vals The values of partials for the reactions
	 * @param i The location on the grid in the depth direction
	 */
	void computeAllPartials(const std::vector<size_t>& startingIdx,
			const std::vector<int>& indices, std::vector<double>& vals,
			int i) const override;
};

} // namespace xolotlCore

#endif
