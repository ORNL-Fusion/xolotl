#ifndef PSISUPERCLUSTER_H
#define PSISUPERCLUSTER_H

// Includes
#include <string>
#include <unordered_map>
#include <Constants.h>
#include "PSICluster.h"
#include "IntegerRange.h"

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
	static std::string buildName(double nHe, double nD, double nT, double nV) {
		std::stringstream nameStream;
		nameStream << "He_" << nHe << "D_" << nD << "T_" << nT << "V_" << nV;
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
		 * The first number represent the moment of A, the second of B
		 * in A + B -> C
		 *
		 * The third number represent which moment we are computing.
		 *
		 * 0 -> l0
		 * 1 -> He
		 * 2 -> D
		 * 3 -> T
		 * 4 -> V
		 */
		double coefs[5][5][5] = { };

		//! The constructor
		ProductionCoefficientBase() {
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
		 * The first number represent the moment of A
		 * in A -> B + C
		 *
		 * The second number represent which moment we are computing.
		 *
		 * 0 -> l0
		 * 1 -> He
		 * 2 -> D
		 * 3 -> T
		 * 4 -> V
		 */
		double coefs[5][5] = { };

		//! The constructor
		SuperClusterDissociationPair(Reaction& _reaction, PSICluster& _first,
				PSICluster& _second) :
				ReactingPairBase(_reaction, _first, _second) {

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

	//! The mean number of atoms in this cluster.
	double numAtom[4] = { };

	//! The total number of clusters gathered in this super cluster.
	int nTot;

	//! The width of the group.
	int sectionWidth[4] = { };

	//! The dispersion in the group.
	double dispersion[4] = { };

	//! The 0th order moment (mean).
	double l0;

	//! The first order moment.
	double l1[4] = { };

	/**
	 * The list of clusters gathered in this.
	 */
	std::set<std::tuple<int, int, int, int> > heVList;

	//! The list of optimized effective reacting pairs.
	ProductionPairMap effReactingList;

	//! The list of optimized effective combining pairs.
	CombiningClusterMap effCombiningList;

	//! The list of optimized effective dissociating pairs.
	DissociationPairMap effDissociatingList;

	//! The list of optimized effective emission pairs.
	DissociationPairMap effEmissionList;

	/**
	 * The first moment flux.
	 */
	double momentFlux[4] = { };

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
	 * @param num The mean number of atoms in this cluster
	 * @param nTot The total number of clusters in this cluster
	 * @param width The width of this super clusteron
	 * @param _network The network
	 * @param registry The performance handler registry
	 */
	PSISuperCluster(double num[4], int nTot, int width[4], IReactionNetwork& _network,
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
	 */
	void resultFrom(ProductionReaction& reaction, int a[4] = defaultInit,
			int b[4] = defaultInit) override;

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
	 */
	void participateIn(ProductionReaction& reaction, int a[4] = defaultInit) override;

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
	 */
	void participateIn(DissociationReaction& reaction, int a[4] = defaultInit,
			int b[4] = defaultInit) override;

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
	 */
	void emitFrom(DissociationReaction& reaction, int a[4] = defaultInit) override;

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
	void setHeVVector(std::set<std::tuple<int, int, int, int> > vec);

	/**
	 * This operation returns the current concentration.
	 *
	 * @param distHe The helium distance in the group
	 * @param distD The deuterium distance in the group
	 * @param distT The tritium distance in the group
	 * @param distV The vacancy distance in the group
	 * @return The concentration of this reactant
	 */
	double getConcentration(double distHe = 0.0, double distD = 0.0,
			double distT = 0.0, double distV = 0.0) const override {
		return l0 + (distHe * l1[0]) + (distD * l1[1]) + (distT * l1[2])
				+ (distV * l1[3]);
	}

	/**
	 * This operation returns the first moment of the given axis.
	 *
	 * @param axis The axis we are intersted in
	 * @return The moment
	 */
	double getMoment(int axis) const override {
		return l1[axis];
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
	 * 0 -> He
	 * 1 -> D
	 * 2 -> T
	 * 3 -> V
	 *
	 * @param atom The number of atoms
	 * @param axis The axis we are intersted in
	 * @return The distance to the mean number of atoms in the group
	 */
	double getDistance(int atom, int axis) const override {
		return (sectionWidth[axis] == 1) ?
				0.0 : 2.0 * (atom - numAtom[axis]) / (sectionWidth[axis] - 1.0);
	}

	/**
	 * This operation returns the factor used for the moments.
	 *
	 * @param atom The number of atoms
	 * @param axis The axis we are intersted in
	 * @return The factor
	 */
	double getFactor(int atom, int axis) const override {
		return (double) (atom - numAtom[axis]) / dispersion[axis];
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
	 * This operation sets the first order moment.
	 *
	 * @param axis The axis we are intersted in
	 * @param mom The moment
	 */
	void setMoment(double mom, int axis) {
		l1[axis] = mom;
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
		momentFlux[0] = 0.0, momentFlux[1] = 0.0, momentFlux[2] = 0.0, momentFlux[3] =
				0.0;

		// Compute the fluxes.
		return getProductionFlux() - getCombinationFlux()
				+ getDissociationFlux() - getEmissionFlux();
	}

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
	 * This operation returns the total change for its first moment.
	 *
	 * @param axis The direction we are interested in
	 * @return The moment flux
	 */
	double getMomentFlux(int axis) const {
		return momentFlux[axis];
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
	void getPartialDerivatives(std::vector<double> & partials) const override;

	/**
	 * This operation computes the partial derivatives due to production
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getProductionPartialDerivatives(std::vector<double> & partials) const
			override;

	/**
	 * This operation computes the partial derivatives due to combination
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getCombinationPartialDerivatives(std::vector<double> & partials) const
			override;

	/**
	 * This operation computes the partial derivatives due to dissociation of
	 * other clusters into this one.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getDissociationPartialDerivatives(std::vector<double> & partials) const
			override;

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getEmissionPartialDerivatives(std::vector<double> & partials) const
			override;

	/**
	 * This operation computes the partial derivatives for the helium moment.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted.
	 * @ param axis The direction
	 */
	void getMomentPartialDerivatives(std::vector<double> & partials, int axis =
			0) const;

	/**
	 * Returns the average number of vacancies.
	 *
	 * @return The average number of vacancies
	 */
	double getNumV() const {
		return numAtom[3];
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
	 * Access bounds on number of given atoms represented by this cluster.
	 *
	 * @ param axis The direction
	 */
	// TODO do we want to make this generic by taking a type parameter?
	const IntegerRange<IReactant::SizeType>& getBounds(int axis) const {
		return bounds[axis];
	}

	/**
	 * Detect if given number of He and V are in this cluster's group.
	 *
	 * @param _nHe number of He of interest.
	 * @param _nV number of V of interest
	 * @return True if _nHe and _nV is contained in our super cluster.
	 */
	bool isIn(IReactant::SizeType nHe, IReactant::SizeType nD,
			IReactant::SizeType nT, IReactant::SizeType nV) const {
		if (!bounds[0].contains(nHe))
			return false;
		if (!bounds[1].contains(nD))
			return false;
		if (!bounds[2].contains(nT))
			return false;
		if (!bounds[3].contains(nV))
			return false;

		return (heVList.find(std::make_tuple(nHe, nD, nT, nV)) != heVList.end());
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

}
/* end namespace xolotlCore */
#endif
