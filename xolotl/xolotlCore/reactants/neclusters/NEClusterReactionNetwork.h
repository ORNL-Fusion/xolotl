#ifndef NE_CLUSTER_REACTION_NETWORK_H
#define NE_CLUSTER_REACTION_NETWORK_H

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
#include "NECluster.h"

namespace xolotlCore {

/**
 *  This class manages the set of reactants and compound reactants (
 *  combinations of normal reactants) for NE clusters. It also manages a
 *  set of properties that describes the total collection.
 *
 *  This class is a very heavyweight class that should not be abused.
 *
 *  Reactants that are added to this network must be added as with shared_ptrs.
 *  Furthermore, reactants that are added to this network have their ids set to
 *  a network specific id. Reactants should not be shared between separate
 *  instances of a NEClusterReactionNetwork.
 */
class NEClusterReactionNetwork: public ReactionNetwork {

private:

	//! The fission rate in nm2/s
	double fissionRate;

	//! The volumetric density of xenon in a bubble in nm-3
	double rho;

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

		return bindingEnergy;
	}

	/**
	 * Add the dissociation connectivity for the reverse reaction if it is allowed.
	 *
	 * @param emittingReactant The reactant that would emit the pair
	 * @param reaction The reaction we want to reverse
	 *
	 */
	void checkForDissociation(IReactant * emittingReactant,
			ProductionReaction& reaction);

	/**
	 * Determine if the reaction is possible given then reactants and product
	 *
	 * @param r1 First reactant.
	 * @param r2 Second reactant.
	 * @param prod Potential product.
	 */
	bool checkOverlap(NECluster& r1, NECluster& r2, NECluster& prod) {
		int width1 = r1.getSectionWidth();
		int size1 = r1.getSize();
		int width2 = r2.getSectionWidth();
		int size2 = r2.getSize();
		int prodWidth = prod.getSectionWidth(), prodSize = prod.getSize();
		int lo1 = ((int) ((double) size1 - (double) width1 / 2.0) + 1), lo2 =
				((int) ((double) size2 - (double) width2 / 2.0) + 1), hi1 =
				((int) ((double) size1 + (double) width1 / 2.0)), hi2 =
				((int) ((double) size2 + (double) width2 / 2.0));
		int prodLo = ((int) ((double) prodSize - (double) prodWidth / 2.0) + 1),
				prodHi = ((int) ((double) prodSize + (double) prodWidth / 2.0));

		int overlap = std::min(prodHi, hi1 + hi2) - std::max(prodLo, lo1 + lo2)
				+ 1;

		if (overlap < 1) return false;
		return true;
	}

public:

	/**
	 * Default constructor, deleted to force construction using parameters.
	 */
	NEClusterReactionNetwork() = delete;

	/**
	 * The Constructor
	 *
	 * @param registry The performance handler registry
	 */
	NEClusterReactionNetwork(
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	NEClusterReactionNetwork(const NEClusterReactionNetwork& other) = delete;

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
	 * This operation returns the number of super reactants in the network.
	 *
	 * @return The number of super reactants in the network
	 */
	int getSuperSize() const override {
		return getAll(ReactantType::NESuper).size();
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
	 * Find the super cluster that contains the original cluster with nHe
	 * helium atoms and nV vacancies.
	 *
	 * @param nXe The number of xenon atoms
	 * @param nD The number of deuterium atoms
	 * @param nT The number of tritium atoms
	 * @param nV The number of vacancies
	 * @return The super cluster representing the cluster with nHe helium
	 * and nV vacancies, or nullptr if no such cluster exists.
	 */
	IReactant * getSuperFromComp(IReactant::SizeType nXe,
			IReactant::SizeType nD, IReactant::SizeType nT,
			IReactant::SizeType nV) const override;

	/**
	 * Get the diagonal fill for the Jacobian, corresponding to the reactions.
	 *
	 * @param diagFill Connectivity map.
	 */
	void getDiagonalFill(SparseFillMap& sfm) override;

	/**
	 * Get the total concentration of atoms contained in the network.
	 *
	 * Here the atoms that are considered are:
	 * 0 xenon
	 *
	 * @param i Index to switch between the different types of atoms
	 * @param minSize The minimum size to take into account
	 * @return The total concentration
	 */
	double getTotalAtomConcentration(int i = 0, int minSize = 0) override;

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
	 * This operation sets the density of xenon in a bubble, needed to compute all the reaction radii
	 * in NE.
	 *
	 * @param density The density
	 */
	void setDensity(double density) override {
		rho = density;
		return;
	}

	/**
	 * This operation returns the density of xenon in a bubble, needed to compute all the reaction radii
	 * in NE.
	 *
	 * @return The density
	 */
	double getDensity() const override {
		return rho;
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

}

#endif
