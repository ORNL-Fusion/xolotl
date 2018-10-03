#ifndef ALLOYCLUSTERREACTIONNETWORK_H
#define ALLOYCLUSTERREACTIONNETWORK_H

// Includes
//#include <xolotlPerf.h>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <ReactionNetwork.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace xolotlCore {

/**
 *  This class manages the set of reactants and compound reactants (
 *  combinations of normal reactants) for Alloy clusters. It also manages a
 *  set of properties that describes the total collection.
 *
 *  This class is a very heavyweight class that should not be abused.
 *
 *  Reactants that are added to this network must be added as with shared_ptrs.
 *  Furthermore, reactants that are added to this network have their ids set to
 *  a network specific id. Reactants should not be shared between separate
 *  instances of a AlloyClusterReactionNetwork.
 */
class AlloyClusterReactionNetwork: public ReactionNetwork {

private:

	/**
	 * Calculate the reaction constant dependent on the
	 * reaction radii and the diffusion coefficients for the
	 * ith and jth clusters, which itself depends on the current
	 * temperature.
	 *
	 * @param reaction The reaction
	 * @param i The location on the grid in the depth direction
	 * @return The rate
	 */
	double calculateReactionRateConstant(const ProductionReaction& reaction,
			int i) const override;

	/**
	 * Calculate the dissociation constant of the first cluster with respect to
	 * the single-species cluster of the same type based on the current clusters
	 * atomic volume, reaction rate constant, and binding energies.
	 *
	 * @param reaction The reaction
	 * @param i The location on the grid in the depth direction
	 * @return The dissociation constant
	 */
	double calculateDissociationConstant(const DissociationReaction& reaction,
			int i) override;
	/**
	 * Calculate the binding energy for the dissociation cluster to emit the single
	 * and second cluster.
	 *
	 * @param reaction The reaction
	 * @return The binding energy corresponding to this dissociation
	 */
	virtual double computeBindingEnergy(
			const DissociationReaction& reaction) const override {
		double bindingEnergy = reaction.first.getFormationEnergy()
				+ reaction.second.getFormationEnergy()
				- reaction.dissociating.getFormationEnergy();

		if (reaction.dissociating.getType() == ReactantType::Void
				|| reaction.dissociating.getType() == ReactantType::VoidSuper) {
			if (reaction.first.getType() == ReactantType::I
					|| reaction.second.getType() == ReactantType::I) {
				double n = reaction.dissociating.getSize();
				bindingEnergy = 3.5
						- 3.45 * (pow(n + 1.0, 2.0 / 3.0) - pow(n, 2.0 / 3.0));
			} else if (reaction.first.getType() == ReactantType::V
					|| reaction.second.getType() == ReactantType::V) {
				double n = reaction.dissociating.getSize();
				bindingEnergy = 1.9
						- 3.45 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
			}
		} else if (reaction.dissociating.getType() == ReactantType::Faulted
				|| reaction.dissociating.getType()
						== ReactantType::FaultedSuper) {
			if (reaction.first.getType() == ReactantType::V
					|| reaction.second.getType() == ReactantType::V) {
				double n = reaction.dissociating.getSize();
				bindingEnergy = 1.9
						- 2.5 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
			}
		} else if (reaction.dissociating.getType() == ReactantType::I) {
			if (reaction.first.getType() == ReactantType::I
					|| reaction.second.getType() == ReactantType::I) {
				double n = reaction.dissociating.getSize();
				bindingEnergy = 3.5
						- 2.5 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
			}
		}

//		std::cout << reaction.dissociating.getName() << " -> " << reaction.first.getName()
//				<< " + " << reaction.second.getName() << std::endl;

		return bindingEnergy;
	}

public:

	/**
	 * Default constructor, deleted to force construction using parameters.
	 */
	AlloyClusterReactionNetwork() = delete;

	/**
	 * The Constructor
	 *
	 * @param registry The performance handler registry
	 */
	AlloyClusterReactionNetwork(
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * The copy constructor.
	 *
	 * @param other
	 */
	AlloyClusterReactionNetwork(const AlloyClusterReactionNetwork &other) = delete;

	/**
	 * To know if the type is vacancy or interstitial
	 *
	 * @param typeName The type of cluster
	 * @return +1 if it is a type of interstitial, -1 if it is a vacancy
	 */
	int typeSwitch(ReactantType const typeName) const;

	/**
	 * Computes the full reaction connectivity matrix for this network.
	 */
	void createReactionConnectivity();

	/**
	 * Add the dissociation connectivity for the reverse reaction if it is allowed.
	 *
	 * @param emittingReactant The reactant that would emit the pair
	 * @param reaction The reaction we want to reverse
	 *
	 */
	void checkDissociationConnectivity(IReactant * emittingReactant,
			std::shared_ptr<ProductionReaction> reaction);

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
	 * This operation returns the number of super reactants in the network.
	 *
	 * @return The number of super reactants in the network
	 */
	int getSuperSize() const override {
		return getAll(ReactantType::VoidSuper).size()
				+ getAll(ReactantType::PerfectSuper).size()
				+ getAll(ReactantType::FrankSuper).size()
				+ getAll(ReactantType::FaultedSuper).size();
	}

	/**
	 * This operation returns the size or number of reactants and momentums in the network.
	 *
	 * @return The number of degrees of freedom
	 */
	virtual int getDOF() const override {
		return size() + getSuperSize() + 1;
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
	 * @param diagFill Connectivity map.
	 */
	void getDiagonalFill(SparseFillMap& sfm) override;

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

}

#endif
