#ifndef NESUPERCLUSTER_H
#define NESUPERCLUSTER_H

// Includes
#include "NECluster.h"
#include <string>
#include <forward_list>

namespace xolotlCore {

/**
 *  A cluster gathering the average properties of many Xe clusters.
 */
class NESuperCluster: public NECluster {

private:
	static std::string buildName(IReactant::SizeType nXe) {
		std::stringstream nameStream;
		nameStream << "Xe_" << nXe;
		return nameStream.str();
	}

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
		NECluster * first;

		/**
		 * The second cluster in the pair
		 */
		NECluster * second;

		/**
		 * The reaction/dissociation pointer to the list
		 */
		Reaction& reaction;

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
		SuperClusterProductionPair(NECluster * firstPtr, NECluster * secondPtr,
				Reaction * _reaction) :
				first(firstPtr), second(secondPtr), reaction(*_reaction), a000(
						0.0), a001(0.0), a010(0.0), a011(0.0), a100(0.0), a101(
						0.0), a110(0.0), a111(0.0) {
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
		NECluster * first;

		/**
		 * The second cluster in the pair
		 */
		NECluster * second;

		/**
		 * The reaction/dissociation pointer to the list
		 */
		Reaction& reaction;

		/**
		 * All the coefficient needed to compute each element
		 */
		double a00;
		double a01;
		double a10;
		double a11;

		//! The constructor
		SuperClusterDissociationPair(NECluster * firstPtr,
				NECluster * secondPtr, Reaction * _reaction) :
				first(firstPtr), second(secondPtr), reaction(*_reaction), a00(
						0.0), a01(0.0), a10(0.0), a11(0.0) {
		}
	};

private:

	//! The mean number of xenon atoms in this cluster.
	double numXe;

	//! The total number of clusters gathered in this super cluster.
	int nTot;

	//! The width in the xenon direction.
	int sectionWidth;

	//! The 0th order moment (mean).
	double l0;

	//! The first order moment in the xenon direction.
	double l1;

	//! The dispersion in the group in the xenon direction.
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
	 * The xenon moment flux.
	 */
	double momentFlux;

public:

	//! The vector of Xe clusters it will replace
	std::vector<NECluster *> xeVector;

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	NESuperCluster() = delete;

	/**
	 * The constructor. All NESuperClusters must be initialized with its
	 * composition.
	 *
	 * @param numXe The mean number of xenon atoms in this cluster
	 * @param nTot The total number of clusters in this cluster
	 * @param width The width of this super cluster in the xenon direction
	 * @param radius The mean radius
	 * @param energy The formation energy
	 * @param _network The network this cluster will belong to.
	 * @param registry The performance handler registry
	 */
	NESuperCluster(double numXe, int nTot, int width, double radius,
			double energy, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	NESuperCluster(NESuperCluster &other) = delete;

	//! Destructor
	~NESuperCluster() {
	}

	/**
	 * Update reactant using other reactants in its network.
	 */
	void updateFromNetwork() override;

	/**
	 * Group the same reactions together and add the reactions to the network lists.
	 */
	void optimizeReactions() override;

	/**
	 * Note that we result from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * \see Reactant.h
	 */
	void resultFrom(ProductionReaction& reaction, double *coef) override;

	/**
	 * Note that we combine with another cluster in a production reaction.
	 * Assumes that the reaction is already in our network.
	 *
	 * \see Reactant.h
	 */
	void participateIn(ProductionReaction& reaction, double *coef) override;

	/**
	 * Note that we combine with another cluster in a dissociation reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * \see Reactant.h
	 */
	void participateIn(DissociationReaction& reaction, double *coef) override;

	/**
	 * Note that we emit from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * \see Reactant.h
	 */
	void emitFrom(DissociationReaction& reaction, double *coef) override;

	/**
	 * This operation returns false.
	 *
	 * @return True if mixed
	 */
	virtual bool isMixed() const override {
		return false;
	}

	/**
	 * This operation returns the total number of clusters it contains.
	 *
	 * @return The total number of clusters
	 */
	int getNTot() const {
		return nTot;
	}

	/**
	 * This operation returns the average number of clusters it contains.
	 *
	 * @return The average number of clusters
	 */
	double getAverage() const {
		return numXe;
	}

	/**
	 * Set the Xe vector
	 */
	void setXeVector(std::vector<NECluster *> vec) {
		xeVector = vec;
	}

	/**
	 * This operation returns the current concentration.
	 *
	 * @param distXe The xenon distance in the group
	 * @return The concentration of this reactant
	 */
	double getConcentration(double distXe) const override {
	    return l0 + (distXe * l1);
    }


	/**
	 * This operation returns the first xenon moment.
	 *
	 * @return The moment
	 */
	double getMoment() const override;

	/**
	 * This operation returns the current total concentration of clusters in the group.

	 * @return The concentration
	 */
	double getTotalConcentration() const;

	/**
	 * This operation returns the current total concentration of xenon in the group.

	 * @return The concentration
	 */
	double getTotalXenonConcentration() const;

	/**
	 * This operation returns the distance to the mean.
	 *
	 * @param xe The number of xenon
	 * @return The distance to the mean number of xenon in the group
	 */
	double getDistance(int xe) const;

	/**
	 * Calculate the dispersion of the group.
	 */
	void computeDispersion();

	/**
	 * Get the dispersion of the group.
	 */
	double getDispersion() const {
		return dispersion;
	}

	/**
	 * This operation sets the zeroth order moment.
	 *
	 * @param mom The moment
	 */
	void setZerothMoment(double mom) {
		l0 = mom;
	}

	/**
	 * This operation sets the first order moment in the xenon direction.
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
	void resetConnectivities() override;

	/**
	 * This operation returns the total flux of this cluster in the
	 * current network.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The total change in flux for this cluster due to all
	 * reactions
	 */
	double getTotalFlux(int i) override;

	/**
	 * This operation returns the total change in this cluster due to
	 * other clusters dissociating into it. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to dissociation of other clusters
	 */
	double getDissociationFlux(int i);

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux(int i);

	/**
	 * This operation returns the total change in this cluster due to
	 * the production of this cluster by other clusters. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to this cluster being produced
	 */
	double getProductionFlux(int i);

	/**
	 * This operation returns the total change in this cluster due to
	 * the combination of this cluster with others. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to this cluster combining with other clusters
	 */
	double getCombinationFlux(int i);

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
	 * @param i The location on the grid in the depth direction
	 *
	 */
	void getPartialDerivatives(std::vector<double> & partials, int i) const
			override;

	/**
	 * This operation computes the partial derivatives due to production
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	void getProductionPartialDerivatives(std::vector<double> & partials,
			int i) const override;

	/**
	 * This operation computes the partial derivatives due to combination
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	void getCombinationPartialDerivatives(std::vector<double> & partials,
			int i) const override;

	/**
	 * This operation computes the partial derivatives due to dissociation of
	 * other clusters into this one.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	void getDissociationPartialDerivatives(std::vector<double> & partials,
			int i) const override;

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	void getEmissionPartialDerivatives(std::vector<double> & partials,
			int i) const override;

	/**
	 * This operation computes the partial derivatives for the xenon moment.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted.
	 */
	void getMomentPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation returns the vector of production reactions in which
	 * this cluster is involved, containing the id of the reactants, and
	 * the a coefs.
	 *
	 * @return The vector of productions
	 */
	virtual std::vector<std::vector<double> > getProdVector() const override;

	/**
	 * This operation returns the vector of combination reactions in which
	 * this cluster is involved, containing the id of the other reactants, and
	 * the a coefs.
	 *
	 * @return The vector of combinations
	 */
	virtual std::vector<std::vector<double> > getCombVector() const override;

	/**
	 * This operation returns the vector of dissociation reactions in which
	 * this cluster is involved, containing the id of the emitting reactants, and
	 * the a coefs.
	 *
	 * @return The vector of dissociations
	 */
	virtual std::vector<std::vector<double> > getDissoVector() const override;

	/**
	 * This operation returns the vector of emission reactions in which
	 * this cluster is involved, containing the a coefs.
	 *
	 * @return The vector of productions
	 */
	virtual std::vector<std::vector<double> > getEmitVector() const override;

	/**
	 * This operation returns the section width.
	 *
	 * @return The width of the section
	 */
	int getSectionWidth() const {
		return sectionWidth;
	}

	/**
	 * Detect if given coordinates are in this cluster's group.
	 *
	 * @param _nXe number of Xe of interest.
	 * @return True if the coordinates are contained in our super cluster.
	 */
	bool isIn(IReactant::SizeType nXe) const {
		return (nXe > numXe - (double) sectionWidth / 2.0
				&& nXe < numXe + (double) sectionWidth / 2.0);
	}

};
//end class NESuperCluster

} /* end namespace xolotlCore */
#endif
