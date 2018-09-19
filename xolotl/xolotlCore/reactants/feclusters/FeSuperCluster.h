#ifndef FESUPERCLUSTER_H
#define FESUPERCLUSTER_H

// Includes
#include <string>
#include <unordered_map>
#include <Constants.h>
#include <MathUtils.h>
#include "FeCluster.h"

namespace xolotlCore {
/**
 *  A cluster gathering the average properties of many HeV clusters.
 */
class FeSuperCluster: public FeCluster {

private:
	static std::string buildName(double nHe, double nV) {
		std::stringstream nameStream;
		nameStream << "He_" << nHe << "V_" << nV;
		return nameStream.str();
	}

protected:

	struct ReactingInfoBase {

		/**
		 * The first cluster in the pair
		 */
		FeCluster& first;

		/**
		 * The reaction/dissociation pointer to the list
		 */
		Reaction& reaction;

		//! The constructor
		ReactingInfoBase(Reaction& _reaction, FeCluster& _first) :
				first(_first), reaction(_reaction) {

		}

		/**
		 * Default and copy constructors, disallowed.
		 */
		ReactingInfoBase() = delete;
		ReactingInfoBase(const ReactingInfoBase& other) = delete;
	};

	struct ReactingPairBase: public ReactingInfoBase {

		/**
		 * The second cluster in the pair
		 */
		FeCluster& second;

		//! The constructor
		ReactingPairBase(Reaction& _reaction, FeCluster& _first,
				FeCluster& _second) :
				ReactingInfoBase(_reaction, _first), second(_second) {

		}

		/**
		 * Default and copy constructors, disallowed.
		 */
		ReactingPairBase() = delete;
		ReactingPairBase(const ReactingPairBase& other) = delete;
	};

	struct ProductionCoefficientBase {

		/**
		 * All the coefficient needed to compute each element
		 * The first number represent the moment of A, the second of B
		 * in A + B -> C
		 *
		 * The third number represent which moment we are computing.
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
		ProductionCoefficientBase() :
				a000(0.0), a001(0.0), a002(0.0), a100(0.0), a101(0.0), a102(
						0.0), a200(0.0), a201(0.0), a202(0.0), a010(0.0), a011(
						0.0), a012(0.0), a020(0.0), a021(0.0), a022(0.0), a110(
						0.0), a111(0.0), a112(0.0), a120(0.0), a121(0.0), a122(
						0.0), a210(0.0), a211(0.0), a212(0.0), a220(0.0), a221(
						0.0), a222(0.0) {
		}

		/**
		 * Copy constructor, disallowed.
		 */
		ProductionCoefficientBase(const ProductionCoefficientBase& other) = delete;
	};

	/**
	 * This is a protected class that is used to implement the flux calculations
	 * for two body production reactions.
	 *
	 * The constants are stored along the clusters taking part in the
	 * reaction or dissociation for faster computation because they only change
	 * when the temperature change. k is computed when setTemperature() is called.
	 */
	struct SuperClusterProductionPair: public ReactingPairBase,
			public ProductionCoefficientBase {

		/**
		 * Nice name for key type in map of key to production pair.
		 */
		using KeyType = ReactantAddrPair;

		//! The constructor
		SuperClusterProductionPair(Reaction& _reaction, FeCluster& _first,
				FeCluster& _second) :
				ReactingPairBase(_reaction, _first, _second), ProductionCoefficientBase() {

		}

		/**
		 * Default and copy constructors, deleted to enforce constructing
		 * using reactants.
		 */
		SuperClusterProductionPair() = delete;
		SuperClusterProductionPair(const SuperClusterProductionPair& other) = delete;
	};

	/**
	 * Concise name for type of map of SuperClusterProductionPairs.
	 */
	using ProductionPairMap = std::unordered_map<SuperClusterProductionPair::KeyType, SuperClusterProductionPair>;

	/**
	 * Info about a cluster we combine with.
	 */
	struct SuperClusterCombiningCluster: public ReactingInfoBase,
			public ProductionCoefficientBase {

		/**
		 * Concise name for type of keys in map of keys to
		 * combining cluster info.
		 */
		using KeyType = IReactant*;

		//! The constructor
		SuperClusterCombiningCluster(Reaction& _reaction, FeCluster& _first) :
				ReactingInfoBase(_reaction, _first), ProductionCoefficientBase() {

		}

		/**
		 * Default and copy construtors, deleted to enforce constructing
		 * using reactants.
		 */
		SuperClusterCombiningCluster() = delete;
		SuperClusterCombiningCluster(const SuperClusterCombiningCluster& other) = delete;
	};

	/**
	 * Concise name for type of map of SuperClusterCombiningClusters.
	 */
	using CombiningClusterMap = std::unordered_map<SuperClusterCombiningCluster::KeyType, SuperClusterCombiningCluster>;

	/**
	 * This is a protected class that is used to implement the flux calculations
	 * for two dissociation reactions.
	 *
	 * The constants are stored along the clusters taking part in the
	 * reaction or dissociation for faster computation because they only change
	 * when the temperature change. k is computed when setTemperature() is called.
	 */
	struct SuperClusterDissociationPair: public ReactingPairBase {

		/**
		 * Concise name for type of key into map of dissociation pairs.
		 */
		using KeyType = ReactantAddrPair;

		/**
		 * All the coefficient needed to compute each element
		 * The first number represent the moment of A
		 * in A -> B + C
		 *
		 * The second number represent which moment we are computing.
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
		SuperClusterDissociationPair(Reaction& _reaction, FeCluster& _first,
				FeCluster& _second) :
				ReactingPairBase(_reaction, _first, _second), a00(0.0), a01(
						0.0), a02(0.0), a10(0.0), a11(0.0), a12(0.0), a20(0.0), a21(
						0.0), a22(0.0) {

		}

		/**
		 * Default and copy constructors, disallowed.
		 */
		SuperClusterDissociationPair() = delete;
		SuperClusterDissociationPair(const SuperClusterDissociationPair& other) = delete;
	};

	/**
	 * Concise name for type of map of SuperClusterDissociationPairs.
	 */
	using DissociationPairMap = std::unordered_map<SuperClusterDissociationPair::KeyType, SuperClusterDissociationPair>;

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

	//! The 0th order moment (mean).
	double l0;

	//! The first order moment in the helium direction.
	double l1He;

	//! The first order moment in the vacancy direction.
	double l1V;

	//! The dispersion in the group in the helium direction.
	double dispersionHe;

	//! The dispersion in the group in the vacancy direction.
	double dispersionV;

	//! The list of optimized effective reacting pairs.
	ProductionPairMap effReactingList;

	//! The list of optimized effective combining pairs.
	CombiningClusterMap effCombiningList;

	//! The list of optimized effective dissociating pairs.
	DissociationPairMap effDissociatingList;

	//! The list of optimized effective emission pairs.
	DissociationPairMap effEmissionList;

	/**
	 * The helium moment flux.
	 */
	double heMomentFlux;

	/**
	 * The vacancy moment flux.
	 */
	double vMomentFlux;

	/**
	 * Output coefficients for a given reaction to the given output stream.
	 *
	 * @param os The output stream on which to write the coefficients.
	 * @param curr Information about our participation in a reaction.
	 */
	void dumpCoefficients(std::ostream& os,
			ProductionCoefficientBase const& curr) const;
	void dumpCoefficients(std::ostream& os,
			SuperClusterDissociationPair const& curr) const;

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	FeSuperCluster() = delete;

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
	FeSuperCluster(double numHe, double numV, int nTot, int heWidth, int vWidth,
			IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	FeSuperCluster(FeSuperCluster &other) = delete;

	//! Destructor
	~FeSuperCluster() {
	}

	/**
	 * Note that we result from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param a Number that can be used by daughter classes.
	 * @param b Number that can be used by daughter classes.
	 */
	void resultFrom(ProductionReaction& reaction, int a[4] = { },
			int b[4] = { }) override;

	/**
	 * Note that we result from the given reaction involving a super cluster.
	 * Assumes the reaction is already in the network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param prInfos Production reaction parameters.
	 */
	void resultFrom(ProductionReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos) override;

	/**
	 * Note that we result from the given reaction involving a super cluster.
	 * Assumes the reaction is already in the network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param product The cluster created by the reaction.
	 */
	void resultFrom(ProductionReaction& reaction, IReactant& product) override;

	/**
	 * Note that we result from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param coef The corresponding coefficient
	 */
	void resultFrom(ProductionReaction& reaction, double *coef) override;

	/**
	 * Note that we combine with another cluster in a production reaction.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param a Number that can be used by daughter classes.
	 */
	void participateIn(ProductionReaction& reaction, int a[4] = { }) override;

	/**
	 * Note that we combine with another cluster in a production reaction
	 * involving a super cluster.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param prInfos Production reaction parameters.
	 */
	void participateIn(ProductionReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos) override;

	/**
	 * Note that we combine with another cluster in a production reaction
	 * involving a super cluster.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param product The cluster created by the reaction.
	 */
	void participateIn(ProductionReaction& reaction, IReactant& product)
			override;

	/**
	 * Note that we combine with another cluster in a production reaction.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param coef Number that can be used by daughter classes.
	 */
	void participateIn(ProductionReaction& reaction, double *coef) override;

	/**
	 * Note that we combine with another cluster in a dissociation reaction.
	 * Assumes the reaction is already inour network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param a Number that can be used by daughter classes.
	 * @param b Number that can be used by daughter classes.
	 */
	void participateIn(DissociationReaction& reaction, int a[4] = { },
			int b[4] = { }) override;

	/**
	 * Note that we combine with another cluster in a dissociation reaction
	 * involving a super cluster.
	 * Assumes the reaction is already inour network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param prInfos Production reaction parameters.
	 */
	void participateIn(DissociationReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos) override;

	/**
	 * Note that we combine with another cluster in a dissociation reaction
	 * involving a super cluster.
	 * Assumes the reaction is already inour network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param disso The dissociating cluster.
	 */
	void participateIn(DissociationReaction& reaction, IReactant& disso)
			override;

	/**
	 * Note that we combine with another cluster in a dissociation reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param coef Number that can be used by daughter classes.
	 */
	void participateIn(DissociationReaction& reaction, double *coef) override;

	/**
	 * Note that we emit from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param a Number that can be used by daughter classes.
	 */
	void emitFrom(DissociationReaction& reaction, int a[4] = { }) override;

	/**
	 * Note that we emit from the given reaction involving a super cluster.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param prInfos Production reaction parameters.
	 */
	void emitFrom(DissociationReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos) override;

	/**
	 * Note that we emit from the given reaction involving a super cluster.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param disso The dissociating cluster.
	 */
	void emitFrom(DissociationReaction& reaction, IReactant& disso) override;

	/**
	 * Note that we emit from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param coef Number that can be used by daughter classes.
	 */
	void emitFrom(DissociationReaction& reaction, double *coef) override;

	/**
	 * This operation returns true to signify that this cluster is a mixture of
	 * He and V.
	 *
	 * @return True if mixed
	 */
	virtual bool isMixed() const override {
		return true;
	}

	/**
	 * Set the HeV vector and compute different parameters
	 */
	void setHeVVector(std::vector<std::pair<int, int> > vec);

	/**
	 * This operation returns the current concentration.
	 *
	 * @param distHe The helium distance in the group
	 * @param distV The vacancy distance in the group
	 * @return The concentration of this reactant
	 */
	double getConcentration(double distHe, double distV) const override {
		return l0 + (distHe * l1He) + (distV * l1V);
	}

	/**
	 * This operation returns the first helium moment.
	 *
	 * @return The moment
	 */
	double getHeMoment() const override {
		return l1He;
	}

	/**
	 * This operation returns the first vacancy moment.
	 *
	 * @return The moment
	 */
	double getVMoment() const override {
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
	double getHeDistance(int he) const override {
		return (sectionHeWidth == 1) ?
				0.0 : 2.0 * (he - numHe) / (sectionHeWidth - 1.0);
	}

	/**
	 * This operation returns the distance to the mean.
	 *
	 * @param he The number of vacancy
	 * @return The distance to the mean number of vacancy in the group
	 */
	double getVDistance(int v) const override {
		return (sectionVWidth == 1) ?
				0.0 : 2.0 * (v - numV) / (sectionVWidth - 1.0);
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
	 * This operation sets the first order moment in the helium direction.
	 *
	 * @param mom The moment
	 */
	void setHeMoment(double mom) {
		l1He = mom;
	}

	/**
	 * This operation sets the first order moment in the vacancy direction.
	 *
	 * @param mom The moment
	 */
	void setVMoment(double mom) {
		l1V = mom;
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
	double getTotalFlux(int i) override {

		// Initialize the fluxes
		heMomentFlux = 0.0;
		vMomentFlux = 0.0;

		// Compute the fluxes.
		return getProductionFlux(i) - getCombinationFlux(i)
				+ getDissociationFlux(i) - getEmissionFlux(i);
	}

	/**
	 * This operation returns the total change in this cluster due to
	 * other clusters dissociating into it. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to dissociation of other clusters
	 */
	double getDissociationFlux(int i) override;

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux(int i) override;

	/**
	 * This operation returns the total change in this cluster due to
	 * the production of this cluster by other clusters. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to this cluster being produced
	 */
	double getProductionFlux(int i) override;

	/**
	 * This operation returns the total change in this cluster due to
	 * the combination of this cluster with others. Compute the contributions to
	 * the moment fluxes at the same time.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to this cluster combining with other clusters
	 */
	double getCombinationFlux(int i) override;

	/**
	 * This operation returns the total change for its helium moment.
	 *
	 * @return The moment flux
	 */
	double getHeMomentFlux() const {
		return heMomentFlux;
	}

	/**
	 * This operation returns the total change for its vacancy moment.
	 *
	 * @return The moment flux
	 */
	double getVMomentFlux() const {
		return vMomentFlux;
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
	 * This operation computes the partial derivatives for the helium moment.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted.
	 */
	void getHeMomentPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives for the vacancy moment.
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
	double getNumV() const {
		return numV;
	}

	/**
	 * Returns the number of clusters contained.
	 *
	 * @return The number of clusters
	 */
	double getNTot() const {
		return nTot;
	}

	/**
	 * Returns the boundaries.
	 *
	 * @return The boundaries
	 */
	Array1D<int, 4> getBounds() const {
		Array1D<int, 4> toReturn;
		toReturn[0] = *(heBounds.begin());
		toReturn[1] = *(heBounds.end());
		toReturn[2] = *(vBounds.begin());
		toReturn[3] = *(vBounds.end());
		return toReturn;
	}

	/**
	 * Detect if given number of He and V are in this cluster's group.
	 *
	 * @param _nHe number of He of interest.
	 * @param _nV number of V of interest
	 * @return True if _nHe and _nV is contained in our super cluster.
	 */
	bool isIn(IReactant::SizeType _nHe, IReactant::SizeType _nV) const {

		return heBounds.contains(_nHe) and vBounds.contains(_nV);
	}

	/**
	 * This operation returns the vector of production reactions in which
	 * this cluster is involved, containing the id of the reactants, the rate, and
	 * the coefs[0][0]
	 *
	 * @return The vector of productions
	 */
	virtual std::vector<std::vector<double> > getProdVector() const override;

	/**
	 * This operation returns the vector of combination reactions in which
	 * this cluster is involved, containing the id of the other reactants, the rate, and
	 * the coefs[0]
	 *
	 * @return The vector of combinations
	 */
	virtual std::vector<std::vector<double> > getCombVector() const override;

	/**
	 * This operation returns the vector of dissociation reactions in which
	 * this cluster is involved, containing the id of the emitting reactants, the rate, and
	 * the coefs[0][0]
	 *
	 * @return The vector of dissociations
	 */
	virtual std::vector<std::vector<double> > getDissoVector() const override;

	/**
	 * This operation returns the vector of emission reactions in which
	 * this cluster is involved, containing the rate, and
	 * the coefs[0][0]
	 *
	 * @return The vector of productions
	 */
	virtual std::vector<std::vector<double> > getEmitVector() const override;

	/**
	 * Tell reactant to output a representation of its reaction coefficients
	 * to the given output stream.
	 *
	 * @param os Output stream on which to output coefficients.
	 */
	virtual void outputCoefficientsTo(std::ostream& os) const override;
};
//end class FeSuperCluster

} /* end namespace xolotlCore */
#endif
