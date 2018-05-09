#ifndef PSISUPERCLUSTER_H
#define PSISUPERCLUSTER_H

// Includes
#include <string>
#include <unordered_map>
#include <cassert>
#include <Constants.h>
#include "PSICluster.h"
#include "ReactionNetwork.h"
#include "IntegerRange.h"
#include "NDArray.h"


// We use std::unordered_map for quick lookup of info about 
// reactions we participate in.
// The C++ standard library defines a std::hash for keys
// that are a single pointer, but not for pairs of pointers,
// so we define our own here.  To improve readability,
// we define a concise name for type of a pair of IReactant pointers
// that we use as keys.

namespace xolotlCore {
/**
 *  A cluster gathering the average properties of many HeV clusters.
 */
class PSISuperCluster: public PSICluster {

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
		PSICluster& first;

		/**
		 * The reaction/dissociation constant associated to this
		 * reaction or dissociation
		 */
		const double& kConstant;

		//! The constructor
		ReactingInfoBase(Reaction& _reaction, PSICluster& _first) :
				first(_first), kConstant(_reaction.kConstant) {

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
		PSICluster& second;

		//! The constructor
		ReactingPairBase(Reaction& _reaction, PSICluster& _first,
				PSICluster& _second) :
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
		 * The first number represent the momentum of A, the second of B
		 * in A + B -> C
		 *
		 * The third number represent which momentum we are computing.
		 *
		 * 0 -> l0
		 * 1 -> He
		 * 2 -> V
		 */
        Array<double, 3, 3, 3> a;

		//! The constructor
		ProductionCoefficientBase() {

            a.Init(0.0);
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
		SuperClusterProductionPair(Reaction& _reaction, PSICluster& _first,
				PSICluster& _second) :
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
		SuperClusterCombiningCluster(Reaction& _reaction, PSICluster& _first) :
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
		 * The first number represent the momentum of A
		 * in A -> B + C
		 *
		 * The second number represent which momentum we are computing.
		 *
		 * 0 -> l0
		 * 1 -> He
		 * 2 -> V
		 */
        Array<double, 3, 3> a;

		//! The constructor
		SuperClusterDissociationPair(Reaction& _reaction, PSICluster& _first,
				PSICluster& _second) :
				ReactingPairBase(_reaction, _first, _second) {

            a.Init(0.0);
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

	/**
	 * Bounds on number of He atoms represented by this cluster.
	 */
	IntegerRange<IReactant::SizeType> heBounds;

	/**
	 * Bounds on number of vacancies represented by this cluster.
	 */
	IntegerRange<IReactant::SizeType> vBounds;

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

	/**
	 * The list of clusters gathered in this.
	 */
	std::set<std::pair<int, int> > heVList;

	//! The list of optimized effective reacting pairs.
	ProductionPairMap effReactingList;

	//! The list of optimized effective combining pairs.
	CombiningClusterMap effCombiningList;

	//! The list of optimized effective dissociating pairs.
	DissociationPairMap effDissociatingList;

	//! The list of optimized effective emission pairs.
	DissociationPairMap effEmissionList;

	/**
	 * The helium momentum flux.
	 */
	double heMomentumFlux;

	/**
	 * The vacancy momentum flux.
	 */
	double vMomentumFlux;

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
	PSISuperCluster() = delete;

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
			int vWidth, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	PSISuperCluster(PSISuperCluster &other) = delete;

	//! Destructor
	~PSISuperCluster() {
	}

	/**
	 * Note that we result from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param a Number that can be used by daughter classes.
	 * @param b Number that can be used by daughter classes.
	 * @param c Number that can be used by daughter classes.
	 * @param d Number that can be used by daughter classes.
	 */
	void resultFrom(ProductionReaction& reaction, int a = 0, int b = 0, int c =
			0, int d = 0) override;

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
	 * Note that we combine with another cluster in a production reaction.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param a Number that can be used by daughter classes.
	 * @param b Number that can be used by daughter classes.
	 */
	void participateIn(ProductionReaction& reaction, int a = 0, int b = 0)
			override;

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
	 * Note that we combine with another cluster in a dissociation reaction.
	 * Assumes the reaction is already inour network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param a Number that can be used by daughter classes.
	 * @param b Number that can be used by daughter classes.
	 * @param c Number that can be used by daughter classes.
	 * @param d Number that can be used by daughter classes.
	 */
	void participateIn(DissociationReaction& reaction, int a = 0, int b = 0,
			int c = 0, int d = 0) override;

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
	 * Note that we emit from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param a Number that can be used by daughter classes.
	 * @param b Number that can be used by daughter classes.
	 * @param c Number that can be used by daughter classes.
	 * @param d Number that can be used by daughter classes.
	 */
	void emitFrom(DissociationReaction& reaction, int a = 0, int b = 0, int c =
			0, int d = 0) override;

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
	void setHeVVector(std::set<std::pair<int, int> > vec);

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
	 * This operation returns the first helium momentum.
	 *
	 * @return The momentum
	 */
	double getHeMomentum() const override {
		return l1He;
	}

	/**
	 * This operation returns the first vacancy momentum.
	 *
	 * @return The momentum
	 */
	double getVMomentum() const override {
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
	 * This operation returns the current concentration for a vacancy number.
	 *
	 * @param v The vacancy number
	 * @return The concentration
	 */
	double getIntegratedVConcentration(int v) const;

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
	void resetConnectivities() override;

	/**
	 * This operation returns the total flux of this cluster in the
	 * current network.
	 *
	 * @return The total change in flux for this cluster due to all
	 * reactions
	 */
	double getTotalFlux() override {
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
	double getHeMomentumFlux() const {
		return heMomentumFlux;
	}

	/**
	 * This operation returns the total change for its vacancy momentum.
	 *
	 * @return The momentum flux
	 */
	double getVMomentumFlux() const {
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
    void computePartialDerivatives(
            double* partials,
            const ReactionNetwork::PartialsIdxMap& partialsIdxMap,
            double* hePartials,
            const ReactionNetwork::PartialsIdxMap& hePartialsIdxMap,
            double* vPartials,
            const ReactionNetwork::PartialsIdxMap& vPartialsIdxMap) const;

	void getPartialDerivatives(std::vector<double> & partials) const override
    {
        assert(false);
    }

	/**
	 * This operation computes the partial derivatives due to production
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void computeProductionPartialDerivatives(
            double* partials,
            const ReactionNetwork::PartialsIdxMap& partialsIdxMap,
            double* hePartials,
            const ReactionNetwork::PartialsIdxMap& hePartialsIdxMap,
            double* vPartials,
            const ReactionNetwork::PartialsIdxMap& vPartialsIdxMap) const;
	void getProductionPartialDerivatives(std::vector<double> & partials) const
			override
    {
        assert(false);
    }

	/**
	 * This operation computes the partial derivatives due to combination
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void computeCombinationPartialDerivatives(
            double* partials,
            const ReactionNetwork::PartialsIdxMap& partialsIdxMap,
            double* hePartials,
            const ReactionNetwork::PartialsIdxMap& hePartialsIdxMap,
            double* vPartials,
            const ReactionNetwork::PartialsIdxMap& vPartialsIdxMap) const;
	void getCombinationPartialDerivatives(std::vector<double> & partials) const
			override
    {
        assert(false);
    }

	/**
	 * This operation computes the partial derivatives due to dissociation of
	 * other clusters into this one.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void computeDissociationPartialDerivatives(
            double* partials,
            const ReactionNetwork::PartialsIdxMap& partialsIdxMap,
            double* hePartials,
            const ReactionNetwork::PartialsIdxMap& hePartialsIdxMap,
            double* vPartials,
            const ReactionNetwork::PartialsIdxMap& vPartialsIdxMap) const;
	void getDissociationPartialDerivatives(std::vector<double> & partials) const
			override
    {
        assert(false);
    }

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void computeEmissionPartialDerivatives(
            double* partials,
            const ReactionNetwork::PartialsIdxMap& partialsIdxMap,
            double* hePartials,
            const ReactionNetwork::PartialsIdxMap& hePartialsIdxMap,
            double* vPartials,
            const ReactionNetwork::PartialsIdxMap& vPartialsIdxMap) const;
	void getEmissionPartialDerivatives(std::vector<double> & partials) const
			override
    {
        assert(false);
    }

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
	 * Access bounds on number of He atoms represented by this cluster.
	 */
	// TODO do we want to make this generic by taking a type parameter?
	const IntegerRange<IReactant::SizeType>& getHeBounds() const {
		return heBounds;
	}

	/**
	 * Access bounds on number of vacancies represented by this cluster.
	 */
	const IntegerRange<IReactant::SizeType>& getVBounds() const {
		return vBounds;
	}

	/**
	 * Detect if given number of He and V are in this cluster's group.
	 *
	 * @param _nHe number of He of interest.
	 * @param _nV number of V of interest
	 * @return True if _nHe and _nV is contained in our super cluster.
	 */
	bool isIn(IReactant::SizeType _nHe, IReactant::SizeType _nV) const {
		if (!heBounds.contains(_nHe))
			return false;
		if (!vBounds.contains(_nV))
			return false;

		return (heVList.find(std::make_pair(_nHe, _nV)) != heVList.end());
	}

	/**
	 * Tell reactant to output a representation of its reaction coefficients
	 * to the given output stream.
	 *
	 * @param os Output stream on which to output coefficients.
	 */
	virtual void outputCoefficientsTo(std::ostream& os) const override;
};
//end class PSISuperCluster

} /* end namespace xolotlCore */
#endif
