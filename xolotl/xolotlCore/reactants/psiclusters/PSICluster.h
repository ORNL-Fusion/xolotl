#ifndef PSICLUSTER_H
#define PSICLUSTER_H

// Includes
#include <Reactant.h>
#include <math.h>
#include <vector>
#include <set>
#include <unordered_map>
#include <sstream>

namespace xolotlCore {

/**
 * The PSICluster class is a Reactant that is specialized to work for
 * simulations of plasma-surface interactions. It provides special routines
 * for calculating the total flux due to production and dissociation and
 * obtaining the cluster size.
 *
 * PSIClusters must always be initialized with a size. If the constructor is
 * passed a size of zero or less, the actual size will be set to 1.
 *
 * The getComposition() operation is implemented by subclasses and will always
 * return a map with the keys He, V, I, HeV or HeI. The operation getTypeName()
 * will always return one of the same values.
 */
class PSICluster: public Reactant {

protected:

	/**
	 * This is a protected class that is used to implement the flux calculations
	 * for two body reactions of the form "first + second -> third".
	 */
	class ReactingPair {
	public:

		/**
		 * The first reacting cluster in the pair
		 */
		PSICluster * first;

		/**
		 * The second reacting cluster in the pair
		 */
		PSICluster * second;

		//! The constructor
		ReactingPair(PSICluster * firstPtr, PSICluster * secondPtr)
		: first(firstPtr), second(secondPtr) {}
	};

	/**
	 * The total size of this cluster including the contributions from all
	 * species.
	 */
	int size;

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
	double diffusionCoefficient;


	/**
	 * The shared pointer to this cluster in the network. It is assigned in
	 * setReactionNetwork().
	 */
	std::shared_ptr<PSICluster> thisSharedPtr;

	/**
	 * The index/id of this cluster in the reaction network - 1. It is used for
	 * indexing arrays (thus the -1).
	 */
	int thisNetworkIndex;

	/**
	 * The binding energies for this cluster with clusters of other species
	 * types. There is one binding energy for each of the other species ordered
	 * by He, V, I and mixed species at indices 0, 1, 2 and 3 respectively.
	 */
	std::vector<double> bindingEnergies;

	/**
	 * A map used to look up binding energy indices quickly.
	 */
	static std::unordered_map<std::string, int> bindingEnergyIndexMap;

	/**
	 * The migration energy for this cluster.
	 */
	double migrationEnergy;

	/**
	 * The reaction radius of this cluster
	 */
	double reactionRadius;

	/**
	 * A vector of ReactingPairs that represent reacting pairs of clusters
	 * that produce this cluster. This vector should be populated early in the
	 * cluster's lifecycle by subclasses. In the standard Xolotl clusters,
	 * this vector is filled in createReactionConnectivity.
	 */
	std::vector<ReactingPair> reactingPairs;

	/**
	 * A vector of clusters that combine with this cluster to produce other
	 * clusters. This vector should be populated early in the cluster's
	 * lifecycle by subclasses. In the standard Xolotl clusters, this vector is
	 * filled in createReactionConnectivity.
	 */
	std::vector<Reactant *> combiningReactants;

	/**
	 * A vector of clusters that dissociate to form this cluster.
	 *
	 * Entries in this vector are added by the dissociating cluster itself, not
	 * the cluster which it creates. The array is filled this way for
	 * mathematical and computational simplicity. For single-species clusters,
	 * this is performed in PSICluster.cpp::dissociateClusters(). Compound
	 * clusters that override dissociateClusters() or don't call it should
	 * make sure that they fill this list.
	 */
	std::vector<PSICluster *> dissociatingClusters;

	/**
	 * A pointer to the cluster of the same type as this one that has
	 * size equal to 1. The reference is set in setReactionNetwork().
	 */
	PSICluster * sameTypeSizeOneCluster;

	/**
	 * Counter for the number of times getDissociationFlux is called.
	 */
	std::shared_ptr<xolotlPerf::IEventCounter> getDissociationFluxCounter;

	/**
	 * Calculate the reaction constant dependent on the
	 * reaction radii and the diffusion coefficients for the
	 * ith and jth clusters, which itself depends on the current
	 * temperature
	 * @param The first cluster interacting
	 * @param The second cluster interacting
	 * @param temperature
	 * @return
	 */
	double calculateReactionRateConstant(const PSICluster & firstcluster,
			const PSICluster & secondcluster, double temperature) const;

	/**
	 * Calculate the dissociation constant of the first cluster with respect to
	 * the single-species cluster of the same type based on the current clusters
	 * atomic volume, reaction rate constant, and binding energies.
	 *
	 * @param One of the clusters that dissociated from the parent
	 * @param The second cluster that dissociated from the parent
	 * @param temperature The current system temperature
	 * @return
	 */
	double calculateDissociationConstant(const PSICluster & firstCluster,
			const PSICluster & secondCluster, double temperature) const;

	/**
	 * Computes a row (or column) of the reaction connectivity matrix
	 * corresponding to this cluster.
	 *
	 * Connections are made between this cluster and any clusters it
	 * affects in combination and production reactions.
	 *
	 * The base-class implementation does nothing and must be overridden
	 * by subclasses.
	 */
	virtual void createReactionConnectivity();

	/**
	 * Computes a row (or column) of the dissociation connectivity matrix
	 * corresponding to this cluster.
	 *
	 * Connections are made between this cluster and any clusters it affects
	 * in a dissociation reaction.
	 *
	 * The base-class implementation handles dissociation for regular clusters
	 * by processing the reaction
	 *
	 * A_x --> A_(x-1) + A
	 *
	 * Compound clusters should implement their own version of this operation.
	 *
	 */
	virtual void createDissociationConnectivity();

	/**
	 * This operation creates the two dissociated clusters from this cluster.
	 * It is called by createDissociationConnectivity to process the reaction
	 * and handle the connectivity.
	 *
	 * @param firstDissociatedCluster The first cluster removed by
	 * dissociation.
	 * @param secondDissociatedCluster The second cluster removed by
	 * dissociation.
	 */
	void dissociateClusters(
			Reactant * firstDissociatedCluster,
			Reactant * secondDissociatedCluster);

	/**
	 * This operation "combines" clusters in the sense that it handles all of
	 * the logic and caching required to correctly process the reaction
	 *
	 * A_x + A_y --> A_(x+y)
	 *
	 * or
	 *
	 * A_x + B_y --> (A_x)(B_y)
	 *
	 * or
	 *
	 * (A_x)(B_y) + B_z --> (A_x)[B_(z+y)]
	 *
	 * for each cluster in the set that interacts with this cluster.
	 *
	 * This operation fills the reaction connectivity array as well as the
	 * array of combining clusters.
	 *
	 * @param clusters The clusters that can combine with this cluster.
	 * @param maxSize The maximum size of the compound produced in the reaction.
	 * @param productName The name of the product produced in the reaction.
	 */
	void combineClusters(std::vector<Reactant *> & clusters,
			int maxSize, std::string productName);

	/**
	 * This operation handles partial replacement reactions of the form
	 *
	 * (A_x)(B_y) + C_z --> (A_x)[B_(y-z)]
	 *
	 * for each compound cluster in the set.
	 *
	 * This operation fills the reaction connectivity array as well as the
	 * array of combining clusters.
	 *
	 * @param clusters The clusters that have part of their B components
	 * replaced. It is assumed that each element of this set represents a
	 * cluster of the form (A_x)(B_y).
	 * @param oldComponentName The name of the component that will be partially
	 * replaced.
	 * @param newComponentName The name of the component that will replace the old
	 * component.
	 */
	void replaceInCompound(std::vector<Reactant *> & clusters,
			std::string oldComponentName, std::string newComponentName);

	/** This operation handles reactions where interstitials fill vacancies,
	 * sometimes referred to vacancy-interstitial annihilation. The reaction
	 * is of the form
	 *
	 * I_a + V_b
	 * --> I_(a-b), if a > b
	 * --> V_(b-a), if a < b
	 * --> 0, if a = b
	 *
	 * It is important to note that I_a + V_b = V_b + I_a.
	 *
	 * The operation assumes "this" is the first cluster in the reaction and
	 * relies on the caller to specify the second cluster name/type.
	 *
	 * This operation fills the reaction connectivity array as well as the
	 * array of combining clusters.
	 *
	 * This operation also I_a and V_b as a reacting pair of the product. It
	 * is simpler and cheaper (O(C)) to do this operation here and quite
	 * computationally difficult by comparison (O(numI*numV)) to do this
	 * operation on the child since it has to search all of the possible
	 * parents.
	 *
	 * @param secondClusterName The name of the first cluster in the reaction,
	 * either "V" or "I" and always the opposite or alternative of
	 * this->getName().
	 * @param clusters The set of clusters of the second type that interact
	 * with this cluster.
	 **/
	void fillVWithI(std::string secondClusterName, std::vector<Reactant *> & clusters);

	/**
	 * This operation computes the partial derivatives due to production
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param temperature The temperature at which the reactions are occurring.
	 */
	virtual void getProductionPartialDerivatives(std::vector<double> & partials, double temperature) const;

	/**
	 * This operation computes the partial derivatives due to combination
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param temperature The temperature at which the reactions are occurring.
	 */
	virtual void getCombinationPartialDerivatives(std::vector<double> & partials, double temperature) const;

	/**
	 * This operation computes the partial derivatives due to dissociation
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param temperature The temperature at which the reactions are occurring.
	 */
	virtual void getDissociationPartialDerivatives(std::vector<double> & partials, double temperature) const;

	/**
	 * This operation prints a forward reaction given the three reactants in
	 * A + B -> C.
	 * @param firstReactant - The first reactant in the reaction, A.
	 * @param secondReactant - The second reactant in the reaction, B.
	 * @param thirdReactant - The third reactant in the reaction, C.
	 */
	void printReaction(const Reactant & firstReactant,
			const Reactant & secondReactant,
			const Reactant & productReactant) const;

	/**
	 * This operation prints a backward reaction given the three reactants in
	 * A -> B + C.
	 * @param firstReactant - The first reactant in the reaction, A.
	 * @param secondReactant - The second reactant in the reaction, B.
	 * @param thirdReactant - The third reactant in the reaction, C.
	 */
	void printDissociation(const Reactant & firstReactant,
			const Reactant & secondReactant,
			const Reactant & productReactant) const;

	/**
	 * This operation signifies that the cluster with cluster Id should be
	 * listed as connected with this cluster through forward reactions.
	 * @param clusterId - The integer id of the cluster that is connected
	 * to this cluster.
	 */
	virtual void setReactionConnectivity(int clusterId);

	/**
	 * This operation returns the connectivity array for this cluster for
	 * forward reactions. An entry with value one means that this cluster
	 * and the cluster with id = index + 1 are connected.
	 * @return The connectivity array for "forward" (non-dissociating)
	 * reactions.
	 */
	virtual std::vector<int> getReactionConnectivity() const;

	/**
	 * This operation returns a set that contains only the entries of the
	 * reaction connectivity array that are non-zero.
	 * @return The set of connected reactants. Each entry in the set is the id
	 * of a connected cluster for forward reactions.
	 */
	std::set<int> getReactionConnectivitySet() const;

	/**
	 * This operation signifies that the cluster with cluster Id should be
	 * listed as connected with this cluster through forward reactions.
	 * @param clusterId - The integer id of the cluster that is connected
	 * to this cluster.
	 */
	virtual void setDissociationConnectivity(int clusterId);

	/**
	 * This operation returns the connectivity array for this cluster for
	 * forward reactions. An entry with value one means that this cluster
	 * and the cluster with id = index + 1 are connected.
	 * @return The connectivity array for "forward" (non-dissociating)
	 * reactions.
	 */
	virtual std::vector<int> getDissociationConnectivity() const;

	/**
	 * This operation returns a set that contains only the entries of the
	 * dissociation connectivity array that are non-zero.
	 * @return The set of connected reactants. Each entry in the set is the id
	 * of a connected cluster for dissociation reactions.
	 */
	const std::set<int> & getDissociationConnectivitySet() const;

	/**
	 * This operation recomputes the diffusion coefficient. It is called
	 * whenever the diffusion factor, migration energy or temperature change.
	 *
	 * @param temp the temperature
	 */
	void recomputeDiffusionCoefficient(double temp);

	/**
	 * This constructor is protected because PSIClusters must always be
	 * initialized with a size.
	 */
	PSICluster(const int clusterSize);

private:

	/**
	 * The row of the reaction connectivity matrix corresponding to
	 * this PSICluster stored as a set.
	 *
	 * If a cluster is involved in a reaction with this PSICluster,
	 * the cluster id is an element of this set.
	 */
	std::set<int> reactionConnectivitySet;

	/**
	 * The row of the dissociation connectivity matrix corresponding to
	 * this PSICluster stored as a set.
	 *
	 * If this PSICluster can dissociate into a particular cluster,
	 * the cluster id is an element of this set.
	 */
	std::set<int> dissociationConnectivitySet;

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	PSICluster();

public:

	/** Constructor
	 * @param clusterSize
	 */
	PSICluster(const int clusterSize,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * The copy constructor
	 * @param other
	 */
	PSICluster(const PSICluster &other);

	/**
	 * The Destructor
	 */
	virtual ~PSICluster();

	/**
	 * This operation returns a cluster that is created using the copy
	 * constructor. If this cluster is actually a subclass of cluster, the
	 * clone will be of the same type and therefore carry all of the members
	 * and virtual functions of the subclass in addition to those of the
	 * cluster. This type of copy is not only handy but, in fact, quite
	 * necessary in those cases where a cluster must be copied but its exact
	 * subclass is unknown and there is no way to make a reasonable assumption
	 * about it.
	 * @return A copy of this cluster.
	 */
	virtual std::shared_ptr<Reactant> clone();

	/**
	 * Sets the collection of other clusters that make up
	 * the reaction network in which this cluster exists.
	 *
	 * @param network The reaction network of which this cluster is a part
	 */
	void setReactionNetwork(
			const std::shared_ptr<ReactionNetwork> reactionNetwork);

	/**
	 * This operation returns the total flux of this cluster in the
	 * current network.
	 * @param temperature The temperature at which to calculate the Diffusion Coefficient
	 * @return The total change in flux for this cluster due to all
	 * reactions
	 */
	virtual double getTotalFlux(const double temperature);

	/**
	 * This operation returns the total change in this cluster due to
	 * dissociation.
	 * @param temperature The temperature at which to calculate the flux
	 * @return The flux due to dissociation.
	 */
	virtual double getDissociationFlux(double temperature) const;

	/**
	 * This operation returns the total change in this cluster due to
	 * the production of this cluster by other clusters.
	 * @param temperature The temperature at which to calculate the flux
	 * @return The flux due to this cluster being produced.
	 */
	virtual double getProductionFlux(double temperature) const;

	/**
	 * This operation returns the total change in this cluster due to
	 * the combination of this cluster with others.
	 * @param temperature The temperature at which to calculate the flux
	 * @return The flux due to this cluster combining with other clusters.
	 */
	virtual double getCombinationFlux(double temperature) const;

	/**
	 * This operation returns the total size of the cluster.
	 * @return The total size of this cluster including the contributions
	 * from all species types.
	 */
	virtual int getSize() const;

	/**
	 * This operation retrieves the binding energy for this cluster.
	 * @return An array of the binding energies of this cluster with clusters
	 * of other types as described above.
	 */
	std::vector<double> getBindingEnergies() const;

	/**
	 * This operation sets the binding energies for this cluster. It expects
	 * the energy vector to be ordered as described above.
	 * @param energies The vector of energies.
	 */
	void setBindingEnergies(std::vector<double> energies);

	/**
	 * This operation retrieves the diffusion factor, D_0, that is used to
	 * calculate the diffusion coefficient for this cluster.
	 * @return The diffusion factor of this cluster
	 */
	double getDiffusionFactor() const;

	/**
	 * This operation sets the diffusion factor, D_0, that is used to calculate
	 * the diffusion coefficient for this cluster.
	 * @param factor The diffusion factor.
	 */
	void setDiffusionFactor(const double factor);

	/**
	 * This operation returns the diffusion coefficient for this cluster and is
	 * calculated from the diffusion factor.
	 * @param temperature The temperature at which to calculate the Diffusion Coefficient
	 * @return The diffusion coefficient.
	 */
	virtual double getDiffusionCoefficient() const;

	/**
	 * This operation sets the migration energy for this cluster.
	 * @param energy The migration energy
	 */
	void setMigrationEnergy(const double energy);

	/**
	 * This operation retrieves the migration energy for this cluster
	 * @return the migration energy
	 */
	double getMigrationEnergy() const;

	/**
	 * This operation returns the reaction radius for the
	 * particular cluster.
	 *
	 * @return
	 */
	virtual double getReactionRadius() const;

	/**
	 * This operation returns a list that represents the connectivity
	 * between this cluster and other clusters in the network.
	 * "Connectivity" indicates whether two clusters interact, via any
	 * mechanism, in an abstract sense (as if they were nodes connected by
	 * an edge on a network graph).
	 *
	 * @return An array of ones and zeros that indicate whether or not this
	 * cluster interacts via any mechanism with another cluster. A "1" at
	 * the i-th entry in this array indicates that the cluster interacts
	 * with the i-th cluster in the ReactionNetwork and a "0" indicates
	 * that it does not.
	 */
	std::vector<int> getConnectivity() const;

	/**
	 * This operation returns the list of partial derivatives of this cluster
	 * with respect to all other clusters in the network. The combined lists
	 * of partial derivatives from all of the clusters in the network can be
	 * used to form, for example, a Jacobian.
	 *
	 * @param the temperature at which the reactions are occurring
	 * @return The partial derivatives for this cluster where index zero
	 * corresponds to the first cluster in the list returned by the
	 * ReactionNetwork::getAll() operation.
	 */
	virtual std::vector<double> getPartialDerivatives(double temperature) const;
	
	/**
	 * This operation works as getPartialDerivatives above, but instead of
	 * returning a vector that it creates it fills a vector that is passed to
	 * it by the caller. This allows the caller to optimize the amount of
	 * memory allocations to just one if they are accessing the partial
	 * derivatives many times.
	 *
	 * @param the temperature at which the reactions are occurring
	 * @param the vector that should be filled with the partial derivatives
	 * for this reactant where index zero corresponds to the first reactant in
	 * the list returned by the ReactionNetwork::getAll() operation. The size of
	 * the vector should be equal to ReactionNetwork::size().
	 *
	 */
	virtual void getPartialDerivatives(double temperature, std::vector<double> & partials) const;

	/**
	 * This operation overrides Reactant's setTemperature operation to
	 * correctly recompute the diffusion coefficient and other
	 * temperature-dependent quantities when the temperature is set.
	 * @param temp
	 */
	virtual void setTemperature(double temp);

	/**
	 * This operation returns true if the cluster is a mixed-species or compound
	 * cluster and false if it is a single species cluster.
	 */
	virtual bool isMixed() const { return false; };

};

} /* end namespace xolotlCore */
#endif
