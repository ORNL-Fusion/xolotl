#ifndef PSICLUSTER_H
#define PSICLUSTER_H

// Includes
#include <Reactant.h>
#include <math.h>
#include <vector>

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
 * return a map with the keys He, V, I, HeV or HeI.
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
		std::shared_ptr<PSICluster> first;
		/**
		 * The second reacting cluster in the pair
		 */
		std::shared_ptr<PSICluster> second;
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
	 * The binding energies for this cluster with clusters of other species
	 * types. There is one binding energy for each of the other species ordered
	 * by He, V, I and mixed species at indices 0, 1, 2 and 3 respectively.
	 */
	std::vector<double> bindingEnergies;

	/**
	 * The migration energy for this cluster.
	 */
	double migrationEnergy;

	/**
	 * The row of the reaction connectivity matrix corresponding to
	 * this PSICluster
	 *
	 * If a cluster is involved in a reaction with this PSICluster,
	 * the element at the cluster's index is 1, otherwise 0.
	 */
	std::vector<int> reactionConnectivity;

	/**
	 * The row of the dissociation connectivity matrix corresponding to
	 * this PSICluster
	 *
	 * If this PSICluster can dissociate into a particular cluster,
	 * the element at the cluster's index is 1, otherwise 0.
	 */
	std::vector<int> dissociationConnectivity;

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
	std::vector<std::shared_ptr<Reactant> > combiningReactants;

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
	std::vector<std::shared_ptr<PSICluster>> dissociatingClusters;

	/**
	 * This operation retrieves the shared_ptr for this cluster from the
	 * network so that the reference count will be maintained.
	 *
	 * Simply creating a shared_ptr with the "this" pointer will break
	 * the reference count, so it is important to have the exact
	 * shared pointer that the network is storing.
	 *
	 * @return The shared_ptr from the network or a null shared_ptr if the
	 * network does not contain this reactant.
	 */
	virtual std::shared_ptr<PSICluster> getThisSharedPtrFromNetwork() const;

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
	 * Return whether or not this PSICluster is a product
	 * of the reaction between clusterI and clusterJ in
	 * this clusters ReactionNetwork. This method should be
	 * specialized by subclasses to indicate whether or not they
	 * are the product of the given reaction.
	 *
	 * @param clusterI
	 * @param clusterJ
	 * @return true if this cluster is a product of i and j
	 */
	virtual bool isProductReactant(const Reactant & clusterI,
			const Reactant & clusterJ);

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
			const std::shared_ptr<Reactant> & firstDissociatedCluster,
			const std::shared_ptr<Reactant> & secondDissociatedCluster);

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
	 * @param compoundName The name of the compound produced in the reaction.
	 */
	void combineClusters(std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>clusters,
			int maxSize, std::string compoundName)
	;

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
	void replaceInCompound(
			std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>clusters,
			std::string oldComponentName, std::string newComponentName);

	/* This operation handles reactions where interstitials fill vacancies,
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
	void fillVWithI(std::string secondClusterName,
			std::shared_ptr<std::vector<std::shared_ptr<Reactant> > > clusters);

private:

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	PSICluster();

public:

	/** Constructor
	 * @param clusterSize
	 */
	PSICluster(const int clusterSize);

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
	virtual double getTotalFlux(double temperature) const;

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
	virtual double getDiffusionCoefficient(double temperature) const;

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
};

} /* end namespace xolotlCore */
#endif
