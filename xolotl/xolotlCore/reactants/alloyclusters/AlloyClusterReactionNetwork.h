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
	 * The map of single-species clusters, indexed by a string representation
	 * of a map that contains the name of the reactant and its size.
	 */
	std::unordered_map<std::string, std::shared_ptr<IReactant> > singleSpeciesMap;

	/**
	 * The map of mixed or compound species clusters, indexed by a
	 * string representation of a map that contains the name of the
	 * constituents of the compound reactant and their sizes.
	 */
	std::unordered_map<std::string, std::shared_ptr<IReactant> > mixedSpeciesMap;

	/**
	 * The map of super species clusters, indexed by a string representation
	 * of a map that contains the name of the constituents of the
	 * compound reactant and their sizes.
	 */
	std::unordered_map<std::string, std::shared_ptr<IReactant> > superSpeciesMap;

	/**
	 * This map stores all of the clusters in the network by type.
	 */
	std::map<std::string,
			std::shared_ptr<std::vector<std::shared_ptr<IReactant> > > > clusterTypeMap;

	/**
	 * Number of void clusters in our network.
	 */
	int numVoidClusters;

	/**
	 * Number of faulted clusters in our network.
	 */
	int numFaultedClusters;

	/**
	 * Number of frank clusters in our network.
	 */
	int numFrankClusters;

	/**
	 * Number of perfect clusters in our network.
	 */
	int numPerfectClusters;

	/**
	 * Maximum size of void clusters in our network.
	 */
	int maxVoidClusterSize;

	/**
	 * Maximum size of faulted clusters in our network.
	 */
	int maxFaultedClusterSize;

	/**
	 * Maximum size of frank clusters in our network.
	 */
	int maxFrankClusterSize;

	/**
	 * Maximum size of perfect clusters in our network.
	 */
	int maxPerfectClusterSize;

	/**
	 * This operation sets the default values of the properties table and names
	 * for this network. It is used on construction and during a copy.
	 */
	void setDefaultPropsAndNames();

	/**
	 * Calculate the reaction constant dependent on the
	 * reaction radii and the diffusion coefficients for the
	 * ith and jth clusters, which itself depends on the current
	 * temperature.
	 *
	 * @param reaction The reaction
	 * @return The rate
	 */
	double calculateReactionRateConstant(ProductionReaction * reaction) const;

	/**
	 * Calculate the dissociation constant of the first cluster with respect to
	 * the single-species cluster of the same type based on the current clusters
	 * atomic volume, reaction rate constant, and binding energies.
	 *
	 * @param reaction The reaction
	 * @return The dissociation constant
	 */
	double calculateDissociationConstant(DissociationReaction * reaction) const;

	/**
	 * The Constructor
	 */
	AlloyClusterReactionNetwork();

public:

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
	AlloyClusterReactionNetwork(const AlloyClusterReactionNetwork &other);

	int typeSwitch(std::string const) const;

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
	 */
	virtual void setTemperature(double temp);

	/**
	 * This operation returns the temperature at which the cluster currently exists.
	 *
	 * @return The temperature.
	 */
	virtual double getTemperature() const;

	/**
	 * This operation returns a reactant with the given type and size if it
	 * exists in the network or null if not.
	 *
	 * @param type the type of the reactant
	 * @param size the size of the reactant
	 * @return A pointer to the reactant
	 */
	IReactant * get(const std::string& type, const int size) const;

	/**
	 * This operation returns a compound reactant with the given type and size
	 * if it exists in the network or null if not.
	 *
	 * @param type the type of the compound reactant
	 * @param sizes an array containing the sizes of each piece of the reactant.
	 * For AlloyClusters, this array must be ordered in size by Xe, V and I. This
	 * array must contain an entry for Xe, V and I, even if only Xe and V or Xe
	 * and I are contained in the mixed-species cluster.
	 * @return A pointer to the compound reactant
	 */
	IReactant * getCompound(const std::string& type,
			const std::vector<int>& sizes) const;

	/**
	 * This operation returns a super reactant with the given type and size
	 * if it exists in the network or null if not.
	 *
	 * @param type The type of the compound reactant
	 * @param size The size of the reactant.
	 * @return A pointer to the super reactant
	 */
	IReactant * getSuper(const std::string& type, const int size) const;

	/**
	 * This operation returns all reactants in the network without regard for
	 * their composition or whether they are compound reactants. The list may
	 * or may not be ordered and the decision is left to implementers.
	 *
	 * @return The list of all of the reactants in the network
	 */
	const std::shared_ptr<std::vector<IReactant *>> & getAll() const;

	/**
	 * This operation returns all reactants in the network with the given name.
	 * The list may or may not be ordered and the decision is left to
	 * implementers.
	 *
	 * @param name The reactant or compound reactant name
	 * @return The list of all of the reactants in the network or null if the
	 * name is invalid.
	 */
	std::vector<IReactant *> getAll(const std::string& name) const;

	/**
	 * This operation adds a reactant or a compound reactant to the network.
	 * Adding a reactant to the network does not set the network as the
	 * reaction network for the reactant. This step must be performed
	 * separately to allow for the scenario where the network is generated
	 * entirely before running.
	 *
	 * This operation sets the id of the reactant to one that is specific
	 * to this network. Do not share reactants across networks! This id is
	 * guaranteed to be between 1 and n, including both, for n reactants in
	 * the network.
	 *
	 * The reactant will not be added to the network if the AlloyCluster does
	 * not recognize it as a type of reactant that it cares about (including
	 * adding null). This operation throws an exception of type std::string
	 * if the reactant is  already in the network.
	 *
	 * @param reactant The reactant that should be added to the network.
	 */
	void add(std::shared_ptr<IReactant> reactant);

	/**
	 * This operation adds a super reactant to the network.
	 * Adding a reactant to the network does not set the network as the
	 * reaction network for the reactant. This step must be performed
	 * separately to allow for the scenario where the network is generated
	 * entirely before running.
	 *
	 * This operation sets the id of the reactant to one that is specific
	 * to this network. Do not share Reactants across networks! This id is
	 * guaranteed to be between 1 and n, including both, for n reactants in
	 * the network.
	 *
	 * The reactant will not be added to the network if the AlloyCluster does
	 * not recognize it as a type of reactant that it cares about (including
	 * adding null). This operation throws an exception of type std::string
	 * if the reactant is  already in the network.
	 *
	 * @param reactant The reactant that should be added to the network.
	 */
	void addSuper(std::shared_ptr<IReactant> reactant);

	/**
	 * This operation removes a group of reactants from the network.
	 *
	 * @param reactants The reactants that should be removed.
	 */
	void removeReactants(const std::vector<IReactant*>& reactants);

	/**
	 * This operation reinitializes the network.
	 *
	 * It computes the cluster Ids and network size from the allReactants vector.
	 */
	void reinitializeNetwork();

	/**
	 * This method redefines the connectivities for each cluster in the
	 * allReactans vector.
	 */
	void reinitializeConnectivities();

	/**
	 * This operation updates the concentrations for all reactants in the
	 * network from an array.
	 *
	 * @param concentrations The array of doubles that will be for the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. Properly aligning the array in memory so that this operation
	 * does not overrun is up to the caller.
	 */
	void updateConcentrationsFromArray(double * concentrations);

	/**
	 * This operation returns the size or number of reactants and momentums in the network.
	 *
	 * @return The number of degrees of freedom
	 */
	virtual int getDOF() {
		return networkSize + numSuperClusters;
	}

	/**
	 * Get the diagonal fill for the Jacobian, corresponding to the reactions.
	 *
	 * @param diagFill The pointer to the vector where the connectivity information is kept
	 */
	void getDiagonalFill(int *diagFill);

	/**
	 * Calculate all the rate constants for the reactions and dissociations of the network.
	 * Need to be called only when the temperature changes.
	 */
	void computeRateConstants();

	/**
	 * Compute the fluxes generated by all the reactions
	 * for all the clusters and their momentums.
	 *
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the fluxes are computed used to find the next solution
	 */
	void computeAllFluxes(double *updatedConcOffset);

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
	virtual void computeAllPartials(double *vals, int *indices, int *size);

	/**
	 * Number of void clusters in our network.
	 */
	int getNumVoidClusters() const {
		return numVoidClusters;
	}

	/**
	 * Number of faulted clusters in our network.
	 */
	int getNumFaultedClusters() const {
		return numFaultedClusters;
	}

	/**
	 * Number of frank clusters in our network.
	 */
	int getNumFrankClusters() const {
		return numFrankClusters;
	}

	/**
	 * Number of perfect clusters in our network.
	 */
	int getNumPerfectClusters() const {
		return numPerfectClusters;
	}

	/**
	 * Maximum size of void clusters in our network.
	 */
	int getMaxVoidClusterSize() const {
		return maxVoidClusterSize;
	}

	/**
	 * Maximum size of faulted clusters in our network.
	 */
	int getMaxFaultedClusterSize() const {
		return maxFaultedClusterSize;
	}

	/**
	 * Maximum size of frank clusters in our network.
	 */
	int getMaxFrankClusterSize() const {
		return maxFrankClusterSize;
	}

	/**
	 * Maximum size of perfect clusters in our network.
	 */
	int getMaxPerfectClusterSize() const {
		return maxPerfectClusterSize;
	}

};

}

#endif
