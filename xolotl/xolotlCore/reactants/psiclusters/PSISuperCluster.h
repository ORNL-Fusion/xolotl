#ifndef PSISUPERCLUSTER_H
#define PSISUPERCLUSTER_H

// Includes
#include "PSICluster.h"
#include <string>
#include <forward_list>

namespace xolotlCore {
/**
 *  A cluster gathering the average properties of many HeV clusters.
 */
class PSISuperCluster: public PSICluster {

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
		PSICluster * first;

		/**
		 * The second cluster in the pair
		 */
		PSICluster * second;

		/**
		 * The reaction/dissociation constant associated to this
		 * reaction or dissociation
		 */
		const double * kConstant;

		/**
		 * All the coefficient needed to compute each element
		 * The first number represent the momentum of A, the second of B
		 * in A + B -> C
		 *
		 * The third number represent which momentum we are computing.
		 *
		 * 0 -> l0
		 * 1 -> He
		 * 2 -> V
		 */
		double a000;
		double a001;
		double a002;
		double a100;
		double a101;
		double a102;
		double a200;
		double a201;
		double a202;
		double a010;
		double a011;
		double a012;
		double a020;
		double a021;
		double a022;
		double a110;
		double a111;
		double a112;
		double a120;
		double a121;
		double a122;
		double a210;
		double a211;
		double a212;
		double a220;
		double a221;
		double a222;

		//! The constructor
		SuperClusterProductionPair(PSICluster * firstPtr,
				PSICluster * secondPtr, Reaction * reaction) :
				first(firstPtr), second(secondPtr), kConstant(
						&(reaction->kConstant)), a000(0.0), a001(0.0), a002(
						0.0), a100(0.0), a101(0.0), a102(0.0), a200(0.0), a201(
						0.0), a202(0.0), a010(0.0), a011(0.0), a012(0.0), a020(
						0.0), a021(0.0), a022(0.0), a110(0.0), a111(0.0), a112(
						0.0), a120(0.0), a121(0.0), a122(0.0), a210(0.0), a211(
						0.0), a212(0.0), a220(0.0), a221(0.0), a222(0.0) {
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
		PSICluster * first;

		/**
		 * The second cluster in the pair
		 */
		PSICluster * second;

		/**
		 * The reaction/dissociation constant associated to this
		 * reaction or dissociation
		 */
		const double * kConstant;

		/**
		 * All the coefficient needed to compute each element
		 * The first number represent the momentum of A
		 * in A -> B + C
		 *
		 * The second number represent which momentum we are computing.
		 *
		 * 0 -> l0
		 * 1 -> He
		 * 2 -> V
		 */
		double a00;
		double a01;
		double a02;
		double a10;
		double a11;
		double a12;
		double a20;
		double a21;
		double a22;

		//! The constructor
		SuperClusterDissociationPair(PSICluster * firstPtr,
				PSICluster * secondPtr, Reaction * reaction) :
				first(firstPtr), second(secondPtr), kConstant(
						&(reaction->kConstant)), a00(0.0), a01(0.0), a02(0.0), a10(
						0.0), a11(0.0), a12(0.0), a20(0.0), a21(0.0), a22(0.0) {
		}
	};

private:

	//! The mean number of helium atoms in this cluster.
	double numHe;

	//! The mean number of atomic vacancies in this cluster.
	double numV;

	//! The total number of clusters gathered in this super cluster.
	int nTot;

	//! The width in the helium direction.
	int sectionHeWidth;

	//! The width in the vacancy direction.
	int sectionVWidth;

	//! The 0th order momentum (mean).
	double l0;

	//! The first order momentum in the helium direction.
	double l1He;

	//! The first order momentum in the vacancy direction.
	double l1V;

	//! The dispersion in the group in the helium direction.
	double dispersionHe;

	//! The dispersion in the group in the vacancy direction.
	double dispersionV;

	//! The map containing all the reacting pairs separated by original composition.
	std::map<std::pair<int, int>, std::vector<ClusterPair> > reactingMap;

	//! The map containing all the combining clusters separated by original composition.
	std::map<std::pair<int, int>, std::vector<CombiningCluster> > combiningMap;

	//! The map containing all the dissociating pairs separated by original composition.
	std::map<std::pair<int, int>, std::vector<ClusterPair> > dissociatingMap;

	//! The map containing all the emission pairs separated by original composition.
	std::map<std::pair<int, int>, std::vector<ClusterPair> > emissionMap;

	//! The list of optimized effective reacting pairs.
	std::forward_list<SuperClusterProductionPair> effReactingList;

	//! The list of optimized effective combining pairs.
	std::forward_list<SuperClusterProductionPair> effCombiningList;

	//! The list of optimized effective dissociating pairs.
	std::forward_list<SuperClusterDissociationPair> effDissociatingList;

	//! The list of optimized effective emission pairs.
	std::forward_list<SuperClusterDissociationPair> effEmissionList;

	/**
	 * The helium momentum flux.
	 */
	double heMomentumFlux;

	/**
	 * The vacancy momentum flux.
	 */
	double vMomentumFlux;

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	PSISuperCluster() :
			PSICluster() {
	}

	/**
	 * Group the same reactions together.
	 */
	void optimizeReactions();

public:

	//! The vector of HeV clusters it will replace
	std::vector<PSICluster *> heVVector;

	/**
	 * The constructor. All SuperClusters must be initialized with its
	 * composition.
	 *
	 * @param numHe The mean number of helium atoms in this cluster
	 * @param numV The mean number of vacancies in this cluster
	 * @param nTot The total number of clusters in this cluster
	 * @param heWidth The width of this super cluster in the helium direction
	 * @param vWidth The width of this super cluster in the vacancy direction
	 * @param radius The mean radius
	 * @param energy The mean formation energy
	 * @param registry The performance handler registry
	 */
	PSISuperCluster(double numHe, double numV, int nTot, int heWidth,
			int vWidth, double radius, double energy,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor.
	 *
	 * @param other the reactant to be copied
	 */
	PSISuperCluster(PSISuperCluster &other);

	//! Destructor
	~PSISuperCluster() {
	}

	/**
	 * This operation returns a Reactant that is created using the copy
	 * constructor of PSISuperCluster.
	 *
	 * @return A copy of this reactant
	 */
	virtual std::shared_ptr<IReactant> clone() {
		return std::make_shared<PSISuperCluster>(*this);
	}

	/**
	 * Sets the collection of other clusters that make up
	 * the reaction network in which this cluster exists.
	 *
	 * @param network The reaction network of which this cluster is a part
	 */
	void setReactionNetwork(
			const std::shared_ptr<IReactionNetwork> reactionNetwork);

	/**
	 * This operation returns true to signify that this cluster is a mixture of
	 * He and V.
	 *
	 * @return True if mixed
	 */
	virtual bool isMixed() const {
		return true;
	}

	/**
	 * Set the HeV vector
	 */
	void setHeVVector(std::vector<PSICluster *> vec) {
		heVVector = vec;
	}

	/**
	 * This operation returns the current concentration.
	 *
	 * @param distHe The helium distance in the group
	 * @param distV The vacancy distance in the group
	 * @return The concentration of this reactant
	 */
	double getConcentration(double distHe, double distV) const {
		return l0 + (distHe * l1He) + (distV * l1V);
	}

	/**
	 * This operation returns the first helium momentum.
	 *
	 * @return The momentum
	 */
	double getHeMomentum() const {
		return l1He;
	}

	/**
	 * This operation returns the first vacancy momentum.
	 *
	 * @return The momentum
	 */
	double getVMomentum() const {
		return l1V;
	}

	/**
	 * This operation returns the current total concentration of clusters in the group.

	 * @return The concentration
	 */
	double getTotalConcentration() const;

	/**
	 * This operation returns the current total concentration of helium in the group.

	 * @return The concentration
	 */
	double getTotalHeliumConcentration() const;

	/**
	 * This operation returns the current total concentration of vacancies in the group.

	 * @return The concentration
	 */
	double getTotalVacancyConcentration() const;

	/**
	 * This operation returns the distance to the mean.
	 *
	 * @param he The number of helium
	 * @return The distance to the mean number of helium in the group
	 */
	double getHeDistance(int he) const {
		return (sectionHeWidth == 1) ?
				0.0 : 2.0 * (he - numHe) / (sectionHeWidth - 1.0);
	}

	/**
	 * This operation returns the distance to the mean.
	 *
	 * @param he The number of vacancy
	 * @return The distance to the mean number of vacancy in the group
	 */
	double getVDistance(int v) const {
		return (sectionVWidth == 1) ?
				0.0 : 2.0 * (v - numV) / (sectionVWidth - 1.0);
	}

	/**
	 * Calculate the dispersion of the group.
	 */
	void computeDispersion();

	/**
	 * This operation sets the zeroth order momentum.
	 *
	 * @param mom The momentum
	 */
	void setZerothMomentum(double mom) {
		l0 = mom;
	}

	/**
	 * This operation sets the first order momentum in the helium direction.
	 *
	 * @param mom The momentum
	 */
	void setHeMomentum(double mom) {
		l1He = mom;
	}

	/**
	 * This operation sets the first order momentum in the vacancy direction.
	 *
	 * @param mom The momentum
	 */
	void setVMomentum(double mom) {
		l1V = mom;
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
	double getTotalFlux() {

		// Initialize the fluxes
		heMomentumFlux = 0.0;
		vMomentumFlux = 0.0;

		// Compute the fluxes.
		return getProductionFlux() - getCombinationFlux()
				+ getDissociationFlux() - getEmissionFlux();
	}

	/**
	 * This operation returns the total change in this cluster due to
	 * other clusters dissociating into it. Compute the contributions to
	 * the momentum fluxes at the same time.
	 *
	 * @return The flux due to dissociation of other clusters
	 */
	double getDissociationFlux();

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation. Compute the contributions to
	 * the momentum fluxes at the same time.
	 *
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux();

	/**
	 * This operation returns the total change in this cluster due to
	 * the production of this cluster by other clusters. Compute the contributions to
	 * the momentum fluxes at the same time.
	 *
	 * @return The flux due to this cluster being produced
	 */
	double getProductionFlux();

	/**
	 * This operation returns the total change in this cluster due to
	 * the combination of this cluster with others. Compute the contributions to
	 * the momentum fluxes at the same time.
	 *
	 * @return The flux due to this cluster combining with other clusters
	 */
	double getCombinationFlux();

	/**
	 * This operation returns the total change for its helium momentum.
	 *
	 * @return The momentum flux
	 */
	double getHeMomentumFlux() {
		return heMomentumFlux;
	}

	/**
	 * This operation returns the total change for its vacancy momentum.
	 *
	 * @return The momentum flux
	 */
	double getVMomentumFlux() {
		return vMomentumFlux;
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
	 * This operation computes the partial derivatives for the helium momentum.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted.
	 */
	void getHeMomentPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives for the vacancy momentum.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted.
	 */
	void getVMomentPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * Returns the average number of vacancies.
	 *
	 * @return The average number of vacancies
	 */
	double getNumV() {
		return numV;
	}

	/**
	 * Returns a vector containing the information about the group's bounderies
	 * in the helium and vacancy directions.
	 *
	 * @return The boundaries
	 */
	std::vector<int> getBoundaries() const {
		std::vector<int> boundaries;
		boundaries.push_back((int) (numHe - (double) sectionHeWidth / 2.0) + 1);
		boundaries.push_back(
				(int) (numHe - (double) sectionHeWidth / 2.0) + sectionHeWidth);
		boundaries.push_back((int) (numV - (double) sectionVWidth / 2.0) + 1);
		boundaries.push_back(
				(int) (numV - (double) sectionVWidth / 2.0) + sectionVWidth);
		return boundaries;
	}

};
//end class PSISuperCluster

} /* end namespace xolotlCore */
#endif
