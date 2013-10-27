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
	 * If a reactant is involved in a reaction with this PSICluster,
	 * the element at the reactant's index is 1, otherwise 0.
	 */
	std::vector<int> reactionConnectivity;

	/**
	 * The row of the dissociation connectivity matrix corresponding to
	 * this PSICluster
	 *
	 * If this PSICluster can dissociate into a particular reactant,
	 * the element at the reactant's index is 1, otherwise 0.
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
	 * A vector of Reactants that combine with this cluster to produce other
	 * clusters. This vector should be populated early in the cluster's
	 * lifecycle by subclasses. In the standard Xolotl clusters, this vector is
	 * filled in createReactionConnectivity.
	 */
	std::vector<std::shared_ptr<Reactant> > combiningReactants;

	/**
	 * Calculate the reaction constant dependent on the
	 * reaction radii and the diffusion coefficients for the
	 * ith and jth reactants, which itself depends on the current
	 * temperature
	 * @param The first reactant interacting
	 * @param The second reactant interacting
	 * @param temperature
	 * @return
	 */
	double calculateReactionRateConstant(const PSICluster & firstReactant,
			const PSICluster & secondReactant, const double temperature) const;

	/**
	 * Calculate the dissociation constant based on the current
	 * reactants atomic volume, reaction rate constant, and binding
	 * energies.
	 *
	 * @param The first reactant interacting
	 * @param The second reactant interacting
	 * @param temperature The current system temperature
	 * @return
	 */
	double calculateDissociationConstant(const PSICluster & firstReactant,
			const PSICluster & secondReactant, const double temperature) const;

	/**
	 * Return whether or not this PSICluster is a product
	 * of the reaction between reactantI and reactantJ in
	 * this Reactants ReactionNetwork. This method should be
	 * specialized by subclasses to indicate whether or not they
	 * are the product of the given reaction.
	 *
	 * @param reactantI
	 * @param reactantJ
	 * @return true if this reactant is a product of i and j
	 */
	virtual bool isProductReactant(const Reactant & reactantI,
			const Reactant & reactantJ);

	/**
	 * Computes a row of the reaction connectivity matrix corresponding to
	 * this reactant.
	 *
	 * If two reactants alone can form a reaction, the element at the position
	 * of the second reactant is 1, otherwise 0.
	 */
	virtual void createReactionConnectivity();

	/**
	 * Computes a row of the dissociation connectivity matrix corresponding to
	 * this reactant.
	 *
	 * If two reactants together can be produced by a single reaction,
	 * the element at the position of the second reactant is 1, otherwise 0.
	 */
	virtual void createDissociationConnectivity();

	/**
	 * This operation "combines" clusters in the sense that it handles all of
	 * the logic and caching required to correctly process the reaction
	 *
	 * A_x + A_y --> A_(x+y)
	 *
	 * A_x + B_y --> (A_x)(B_y)
	 *
	 * or
	 *
	 * (A_x)(B_y) + B_z --> (A_x)[B_(z+y)]
	 *
	 * for each reactant in the set that interacts with this reactant.
	 *
	 * @param reactants The reactants that can combine with this reactant.
	 * @param maxSize The maximum size of the compound produced in the reaction.
	 * @param compoundName The name of the compound produced in the reaction.
	 */
	void combineClusters(std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>reactants,
			int maxSize, std::string compoundName);

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
	 * This operation returns a Reactant that is created using the copy
	 * constructor. If this Reactant is actually a subclass of Reactant, the
	 * clone will be of the same type and therefore carry all of the members
	 * and virtual functions of the subclass in addition to those of the
	 * Reactant. This type of copy is not only handy but, in fact, quite
	 * necessary in those cases where a Reactant must be copied but its exact
	 * subclass is unknown and there is no way to make a reasonable assumption
	 * about it.
	 * @return A copy of this reactant.
	 */
	virtual std::shared_ptr<Reactant> clone();

	/**
	 * Sets the collection of other reactants that make up
	 * the reaction network in which this reactant exists.
	 *
	 * @param network The reaction network of which this reactant is a part
	 */
	void setReactionNetwork(
			const std::shared_ptr<ReactionNetwork> reactionNetwork);

	/**
	 * This operation returns the total flux of this reactant in the
	 * current network.
	 * @param temperature The temperature at which to calculate the Diffusion Coefficient
	 * @return The total change in flux for this reactant due to all
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
	 * This virtual method is for subclasses to specialize
	 * to return their representative cluster map, which is a mapping
	 * of which species exist in the cluster to the integer number
	 * of each species.
	 *
	 * @return
	 */
	virtual std::map<std::string, int> getClusterMap();

	/**
	 * This virtual method is for subclasses to specialize
	 * and should return the reaction radius for the
	 * particular PSICluster subclass.
	 *
	 * @return
	 */
	virtual const double getReactionRadius() const;

	/**
	 * This operation returns a list that represents the connectivity
	 * between this Reactant and other Reactants in the network.
	 * "Connectivity" indicates whether two Reactants interact, via any
	 * mechanism, in an abstract sense (as if they were nodes connected by
	 * an edge on a network graph).
	 *
	 * @return An array of ones and zeros that indicate whether or not this
	 * Reactant interacts via any mechanism with another Reactant. A "1" at
	 * the i-th entry in this array indicates that the Reactant interacts
	 * with the i-th Reactant in the ReactionNetwork and a "0" indicates
	 * that it does not.
	 */
	std::vector<int> getConnectivity() const;

	/**
	 * This operation returns the list of partial derivatives of this Reactant
	 * with respect to all other reactants in the network. The combined lists
	 * of partial derivatives from all of the reactants in the network can be
	 * used to form, for example, a Jacobian.
	 *
	 * @param the temperature at which the reactions are occurring
	 * @return The partial derivatives for this reactant where index zero
	 * corresponds to the first reactant in the list returned by the
	 * ReactionNetwork::getAll() operation.
	 */
	std::vector<double> getPartialDerivatives(double temperature) const;
};

} /* end namespace xolotlCore */
#endif
