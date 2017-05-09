#ifndef ALLOYSUPERCLUSTER_H
#define ALLOYSUPERCLUSTER_H

// Includes
#include "AlloyCluster.h"
#include <string>
#include <forward_list>

namespace xolotlCore {

/**
 *  A cluster gathering the average properties of many atom clusters.
 *  Atom here can represent different types.
 */
class AlloySuperCluster: public AlloyCluster {

protected:

	/**
	 * This is a protected class that is used to implement the flux calculations
	 * for two body production reactions.
	 *
	 * The constants are stored along the clusters taking part in the
	 * reaction or dissociation for faster computation because they only change
	 * when the temperature change. k is computed when setTemperature() is called.
	 */
	class SuperClusterProductionPair {
	public:

		/**
		 * The first cluster in the pair
		 */
		AlloyCluster * first;

		/**
		 * The second cluster in the pair
		 */
		AlloyCluster * second;

		/**
		 * The reaction/dissociation constant associated to this
		 * reaction or dissociation
		 */
		const double * kConstant;

		/**
		 * All the coefficient needed to compute each element
		 */
		double a000;
		double a001;
		double a010;
		double a011;
		double a100;
		double a101;
		double a110;
		double a111;

		//! The constructor
		SuperClusterProductionPair(AlloyCluster * firstPtr, AlloyCluster * secondPtr,
				Reaction * reaction) :
				first(firstPtr), second(secondPtr), kConstant(
						&(reaction->kConstant)), a000(0.0), a001(0.0), a010(
						0.0), a011(0.0), a100(0.0), a101(0.0), a110(0.0), a111(
						0.0) {
		}
	};

	/**
	 * This is a protected class that is used to implement the flux calculations
	 * for two dissociation reactions.
	 *
	 * The constants are stored along the clusters taking part in the
	 * reaction or dissociation for faster computation because they only change
	 * when the temperature change. k is computed when setTemperature() is called.
	 */
	class SuperClusterDissociationPair {
	public:

		/**
		 * The first cluster in the pair
		 */
		AlloyCluster * first;

		/**
		 * The second cluster in the pair
		 */
		AlloyCluster * second;

		/**
		 * The reaction/dissociation constant associated to this
		 * reaction or dissociation
		 */
		const double * kConstant;

		/**
		 * All the coefficient needed to compute each element
		 */
		double a00;
		double a01;
		double a10;
		double a11;

		//! The constructor
		SuperClusterDissociationPair(AlloyCluster * firstPtr,
				AlloyCluster * secondPtr, Reaction * reaction) :
				first(firstPtr), second(secondPtr), kConstant(
						&(reaction->kConstant)), a00(0.0), a01(0.0), a10(0.0), a11(
						0.0) {
		}
	};

private:

	//! The mean number of atoms in this cluster.
	double numAtom;

	//! The total number of clusters gathered in this super cluster.
	int nTot;

	//! The width.
	int sectionWidth;

	//! The 0th order moment (mean).
	double l0;

	//! The first order moment.
	double l1;

	//! The dispersion in the group.
	double dispersion;

	//! The map containing all the reacting pairs separated by original composition.
	std::map<int, std::vector<ClusterPair> > reactingMap;

	//! The map containing all the combining clusters separated by original composition.
	std::map<int, std::vector<CombiningCluster> > combiningMap;

	//! The map containing all the dissociating pairs separated by original composition.
	std::map<int, std::vector<ClusterPair> > dissociatingMap;

	//! The map containing all the emission pairs separated by original composition.
	std::map<int, std::vector<ClusterPair> > emissionMap;

	//! The list of optimized effective reacting pairs.
	std::forward_list<SuperClusterProductionPair> effReactingList;

	//! The list of optimized effective combining pairs.
	std::forward_list<SuperClusterProductionPair> effCombiningList;

	//! The list of optimized effective dissociating pairs.
	std::forward_list<SuperClusterDissociationPair> effDissociatingList;

	//! The list of optimized effective emission pairs.
	std::forward_list<SuperClusterDissociationPair> effEmissionList;

	/**
	 * The moment flux.
	 */
	double momentFlux;

	/**
	 * The default constructor is private because AlloyClusters must always be
	 * initialized with a size.
	 */
	AlloySuperCluster() :
			AlloyCluster() {
	}

public:

	//! The vector of atom clusters it will replace
	std::vector<AlloyCluster *> atomVector;

	/**
	 * The constructor. All AlloySuperClusters must be initialized with its
	 * composition.
	 *
	 * @param numAtom The mean number of atoms in this cluster
	 * @param nTot The total number of clusters in this cluster
	 * @param width The width of this super cluster
	 * @param radius The mean radius
	 * @param energy The formation energy
	 * @param type The type of cluster
	 * @param registry The performance handler registry
	 */
	AlloySuperCluster(double numAtom, int nTot, int width, double radius,
			double energy, std::string type,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor.
	 *
	 * @param other the reactant to be copied
	 */
	AlloySuperCluster(AlloySuperCluster &other);

	//! Destructor
	~AlloySuperCluster() {
	}

	/**
	 * This operation returns a Reactant that is created using the copy
	 * constructor of AlloySuperCluster.
	 *
	 * @return A copy of this reactant
	 */
	virtual std::shared_ptr<IReactant> clone();

	/**
	 * Sets the collection of other clusters that make up
	 * the reaction network in which this cluster exists.
	 *
	 * @param network The reaction network of which this cluster is a part
	 */
	void setReactionNetwork(
			const std::shared_ptr<IReactionNetwork> reactionNetwork);

	/**
	 * Group the same reactions together and add the reactions to the network lists.
	 */
	void optimizeReactions();

	/**
	 * This operation returns false.
	 *
	 * @return True if mixed
	 */
	virtual bool isMixed() const {
		return false;
	}

	/**
	 * Set the atom vector
	 */
	void setAtomVector(std::vector<AlloyCluster *> vec) {
		atomVector = vec;
	}

	/**
	 * This operation returns the current concentration.
	 *
	 * @param distAtom The distance in the group
	 * @param distB Unused here
	 * @return The concentration of this reactant
	 */
	double getConcentration(double distAtom, double distB = 0.0) const;

	/**
	 * This operation returns the first moment.
	 *
	 * @return The moment
	 */
	double getMoment() const;

	/**
	 * This operation returns the current total concentration of clusters in the group.

	 * @return The concentration
	 */
	double getTotalConcentration() const;

	/**
	 * This operation returns the current total concentration of atoms in the group.

	 * @return The concentration
	 */
	double getTotalAtomConcentration() const;

	/**
	 * This operation returns the distance to the mean.
	 *
	 * @param atom The number of atom
	 * @return The distance to the mean number of in the group
	 */
	double getDistance(int atom) const;

	/**
	 * Calculate the dispersion of the group.
	 */
	void computeDispersion();

	/**
	 * This operation sets the zeroth order moment.
	 *
	 * @param mom The moment
	 */
	void setZerothMoment(double mom) {
		l0 = mom;
	}

	/**
	 * This operation sets the first order moment.
	 *
	 * @param mom The moment
	 */
	void setMoment(double mom) {
		l1 = mom;
	}

	/**
	 * This operation reset the connectivity sets based on the information
	 * in the production and dissociation vectors.
	 */
	void resetConnectivities();

	/**
	 * This operation returns the total flux of this cluster in the
	 * current network.
	 *
	 * @return The total change in flux for this cluster due to all
	 * reactions
	 */
	double getTotalFlux();

	/**
	 * This operation returns the total change in this cluster due to
	 * other clusters dissociating into it. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @return The flux due to dissociation of other clusters
	 */
	double getDissociationFlux();

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux();

	/**
	 * This operation returns the total change in this cluster due to
	 * the production of this cluster by other clusters. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @return The flux due to this cluster being produced
	 */
	double getProductionFlux();

	/**
	 * This operation returns the total change in this cluster due to
	 * the combination of this cluster with others. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @return The flux due to this cluster combining with other clusters
	 */
	double getCombinationFlux();

	/**
	 * This operation returns the total change for its moment.
	 *
	 * @return The moment flux
	 */
	double getMomentFlux() {
		return momentFlux;
	}

	/**
	 * This operation works as getPartialDerivatives above, but instead of
	 * returning a vector that it creates it fills a vector that is passed to
	 * it by the caller. This allows the caller to optimize the amount of
	 * memory allocations to just one if they are accessing the partial
	 * derivatives many times.
	 *
	 * @param the vector that should be filled with the partial derivatives
	 * for this reactant where index zero corresponds to the first reactant in
	 * the list returned by the ReactionNetwork::getAll() operation. The size of
	 * the vector should be equal to ReactionNetwork::size().
	 *
	 */
	void getPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives due to production
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getProductionPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives due to combination
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getCombinationPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives due to dissociation of
	 * other clusters into this one.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getDissociationPartialDerivatives(
			std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getEmissionPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives for the moment.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted.
	 */
	void getMomentPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation returns the section width.
	 *
	 * @return The width of the section
	 */
	int getSectionWidth() const {
		return sectionWidth;
	}

};
//end class AlloySuperCluster

} /* end namespace xolotlCore */
#endif
