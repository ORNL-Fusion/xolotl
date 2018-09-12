#ifndef FECLUSTER_H
#define FECLUSTER_H

// Includes
#include <Reactant.h>
#include "IntegerRange.h"

namespace xolotlPerf {
class ITimer;
}

namespace xolotlCore {

/**
 * The FeCluster class is a Reactant that is specialized to work for
 * simulations of irradiated iron. It provides special routines
 * for calculating the total flux due to production and dissociation and
 * obtaining the cluster size.
 *
 * FeClusters must always be initialized with a size. If the constructor is
 * passed a size of zero or less, the actual size will be set to 1.
 *
 * The getComposition() operation is implemented by subclasses and will always
 * return a map with the keys He, V, I, HeV. The operation getTypeName()
 * will always return one of the same values.
 *
 * As a rule, it is possible to access directly some of the private members of
 * this class (id, concentration, reactionRadius, diffusionCoefficient, size,
 * type) instead of using the "get" functions for performance reasons. In
 * order to change these values the "set" functions must still be used.
 */
class FeCluster: public Reactant {

protected:

	/**
	 * This is a protected class that is used to implement the flux calculations
	 * for two body reactions or dissociation.
	 *
	 * The constant k+ or k- is stored along the clusters taking part in the
	 * reaction or dissociation for faster computation because they only change
	 * when the temperature change. k is computed when setTemperature() is called.
	 */
	class ClusterPair {
	public:

		/**
		 * The first cluster in the pair
		 */
		FeCluster& first;

		/**
		 * The second cluster in the pair
		 */
		FeCluster& second;

		/**
		 * The reaction/dissociation pointer to the list
		 */
		Reaction& reaction;

		/**
		 * All the coefficient needed to compute each element
		 * The first number represent the moment of A, the second of B
		 * in A + B -> C
		 *
		 * 0 -> l0
		 * 1 -> He
		 * 2 -> V
		 */
		double a00;
		double a10;
		double a20;
		double a01;
		double a02;
		double a11;
		double a12;
		double a21;
		double a22;

		//! The constructor
		ClusterPair(Reaction& _reaction, FeCluster& _first, FeCluster& _second) :
				first(_first), second(_second), reaction(_reaction), a00(0.0), a10(
						0.0), a20(0.0), a01(0.0), a02(0.0), a11(0.0), a12(0.0), a21(
						0.0), a22(0.0) {
		}

		/**
		 * Default and copy constructors, disallowed.
		 */
		ClusterPair() = delete;

		// NB: if FeCluster keeps these in a std::vector,
		// copy ctor is needed.  Implicit definition is fine.
		// CombiningCluster(const CombiningCluster& other) = delete;
		// ClusterPair(const ClusterPair& other) = delete;
	};

	/**
	 * This is a protected class that is used to implement the flux calculations
	 * for combinations.
	 *
	 * The constant k+ is stored along the cluster that combines with this cluster
	 * for faster computation because they only change when the temperature change.
	 * k+ is computed when setTemperature() is called.
	 */
	struct CombiningCluster {

		/**
		 * The combining cluster
		 */
		FeCluster& combining;

		/**
		 * The reaction pointer to the list
		 */
		Reaction& reaction;

		/**
		 * All the coefficient needed to compute each element
		 * The first number represent the moment of A
		 * in A + this -> C
		 *
		 * 0 -> l0
		 * 1 -> He
		 * 2 -> V
		 */
		double a0;
		double a1;
		double a2;

		//! The constructor
		CombiningCluster(Reaction& _reaction, FeCluster& _comb) :
				combining(_comb), reaction(_reaction), a0(0.0), a1(0.0), a2(0.0) {
		}

		/**
		 * Default constructor, disallowed to prohibit building without args.
		 */
		CombiningCluster() = delete;

		// NB: if FeCluster keeps these in a std::vector,
		// copy ctor is needed.  Implicit definition is fine.
		// CombiningCluster(const CombiningCluster& other) = delete;
	};

	/**
	 * Bounds on number of He atoms represented by this cluster.
	 */
	IntegerRange<IReactant::SizeType> heBounds;

	/**
	 * Bounds on number of vacancies represented by this cluster.
	 */
	IntegerRange<IReactant::SizeType> vBounds;

	/**

	 * This operation returns a set that contains only the entries of the
	 * reaction connectivity array that are non-zero.
	 *
	 * @return The set of connected reactants. Each entry in the set is the id
	 * of a connected cluster for forward reactions.
	 */
	std::set<int> getReactionConnectivitySet() const;

	/**
	 * This operation returns a set that contains only the entries of the
	 * dissociation connectivity array that are non-zero.
	 *
	 * @return The set of connected reactants. Each entry in the set is the id
	 * of a connected cluster for dissociation reactions
	 */
	const std::set<int> & getDissociationConnectivitySet() const;

	/**
	 * Output coefficients for a given reaction to the given output stream.
	 *
	 * @param os The output stream on which to write the coefficients.
	 * @param curr Information about our participation in a reaction.
	 */
	void dumpCoefficients(std::ostream& os, ClusterPair const& curr) const;
	void dumpCoefficients(std::ostream& os, CombiningCluster const& curr) const;

public:

	/**
	 * A vector of ClusterPairs that represents reacting pairs of clusters
	 * that produce this cluster. This vector should be populated early in the
	 * cluster's lifecycle by subclasses. In the standard Xolotl clusters,
	 * this vector is filled in createReactionConnectivity.
	 */
	std::vector<ClusterPair> reactingPairs;

	/**
	 * A vector of clusters that combine with this cluster to produce other
	 * clusters. This vector should be populated early in the cluster's
	 * lifecycle by subclasses. In the standard Xolotl clusters, this vector is
	 * filled in createReactionConnectivity.
	 */
	std::vector<CombiningCluster> combiningReactants;

	/**
	 * A vector of pairs of clusters: the first one is the one dissociation into
	 * this cluster, the second one is the one that is emitted at the same time
	 * during the dissociation. This vector should be populated early in the
	 * cluster's lifecycle by subclasses. In the standard Xolotl clusters, this
	 * vector is filled in dissociateCluster that is called by
	 * createDissociationConnectivity.
	 */
	std::vector<ClusterPair> dissociatingPairs;

	/**
	 * A vector of ClusterPairs that represent pairs of clusters that are emitted
	 * from the dissociation of this cluster. This vector should be populated early
	 * in the cluster's lifecycle by subclasses. In the standard Xolotl clusters,
	 * this vector is filled in emitClusters that is called by
	 * createDissociationConnectivity.
	 */
	std::vector<ClusterPair> emissionPairs;

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	FeCluster() = delete;

	/**
	 * Construct a FeCluster.
	 *
	 * @param registry The performance handler registry
	 */
	FeCluster(IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry,
			const std::string& _name = "FeCluster") :
			Reactant(_network, registry, _name), heBounds(0, 0), vBounds(0, 0) {

	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	FeCluster(FeCluster &other) = delete;

	/**
	 * The destructor
	 */
	virtual ~FeCluster() {
	}

	/**
	 * Update reactant using other reactants in its network.
	 */
	virtual void updateFromNetwork() override;

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
	 * @param product The product cluster.
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
	 * @param product The product cluster.
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
	 * This operation returns the connectivity array for this cluster for
	 * forward reactions. An entry with value one means that this cluster
	 * and the cluster with id = index + 1 are connected.
	 * 
	 * @return The connectivity array for "forward" (non-dissociating)
	 * reactions
	 */
	virtual std::vector<int> getReactionConnectivity() const;

	/**
	 * This operation returns the connectivity array for this cluster for
	 * forward reactions. An entry with value one means that this cluster
	 * and the cluster with id = index + 1 are connected.
	 * 
	 * @return The connectivity array for "forward" (non-dissociating)
	 * reactions
	 */
	virtual std::vector<int> getDissociationConnectivity() const;

	/**
	 * This operation returns the current concentration.
	 *
	 * @param distHe The helium distance in the group
	 * @param distV The vacancy distance in the group
	 * @return The concentration of this reactant
	 */
	virtual double getConcentration(double distHe, double distV) const {
        // TODO should this version ever be called?  It ignores
        // its arguments.
		return concentration;
	}

    /**
     * Obtain current concentration.
     *
     * We hid the base class' zero-argument getConcentration when
     * we defined the two-argument version here.  This tells the compiler
     * to allow us to use the base class' zero-argument version also
     * without having to explicitly specify the class scope when calling it.
     *
     * @return The concentration of the reactant.
     */
    using Reactant::getConcentration;

	/**
	 * This operation returns the first helium moment.
	 *
	 * @return The moment
	 */
	virtual double getHeMoment() const {
		return 0.0;
	}

	/**
	 * This operation returns the first vacancy moment.
	 *
	 * @return The moment
	 */
	virtual double getVMoment() const {
		return 0.0;
	}
	/**
	 * This operation returns the distance to the mean.
	 *
	 * @param he The number of helium
	 * @return The distance to the mean number of helium in the group
	 */
	virtual double getHeDistance(int he) const {
		return 0.0;
	}

	/**
	 * This operation returns the distance to the mean.
	 *
	 * @param he The number of vacancy
	 * @return The distance to the mean number of vacancy in the group
	 */
	virtual double getVDistance(int v) const {
		return 0.0;
	}

	/**
	 * This operation returns the total flux of this cluster in the
	 * current network.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The total change in flux for this cluster due to all
	 * reactions
	 */
	virtual double getTotalFlux(int i) override {
		return getProductionFlux(i) - getCombinationFlux(i)
				+ getDissociationFlux(i) - getEmissionFlux(i);
	}

	/**
	 * This operation returns the total change in this cluster due to
	 * other clusters dissociating into it.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to dissociation of other clusters
	 */
	virtual double getDissociationFlux(int i) const;

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to its dissociation
	 */
	virtual double getEmissionFlux(int i) const;

	/**
	 * This operation returns the total change in this cluster due to
	 * the production of this cluster by other clusters.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to this cluster being produced
	 */
	virtual double getProductionFlux(int i) const;

	/**
	 * This operation returns the total change in this cluster due to
	 * the combination of this cluster with others.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The flux due to this cluster combining with other clusters
	 */
	virtual double getCombinationFlux(int i) const;

	/**
	 * This operation returns the list of partial derivatives of this cluster
	 * with respect to all other clusters in the network. The combined lists
	 * of partial derivatives from all of the clusters in the network can be
	 * used to form, for example, a Jacobian.
	 *
	 * @param i The location on the grid in the depth direction
	 * @return The partial derivatives for this cluster where index zero
	 * corresponds to the first cluster in the list returned by the
	 * ReactionNetwork::getAll() operation.
	 */
	virtual std::vector<double> getPartialDerivatives(int i) const override;

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
	virtual void getPartialDerivatives(std::vector<double> & partials,
			int i) const override;

	/**
	 * This operation computes the partial derivatives due to production
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	virtual void getProductionPartialDerivatives(std::vector<double> & partials,
			int i) const;

	/**
	 * This operation computes the partial derivatives due to combination
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	virtual void getCombinationPartialDerivatives(
			std::vector<double> & partials, int i) const;

	/**
	 * This operation computes the partial derivatives due to dissociation of
	 * other clusters into this one.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	virtual void getDissociationPartialDerivatives(
			std::vector<double> & partials, int i) const;

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	virtual void getEmissionPartialDerivatives(std::vector<double> & partials,
			int i) const;

	/**
	 * This operation reset the connectivity sets based on the information
	 * in the effective production and dissociation vectors.
	 */
	void resetConnectivities() override;

	/**
	 * This operation returns the sum of combination rate and emission rate
	 * (where this cluster is on the left side of the reaction) for this
	 * particular cluster.
	 * This is used to computed the desorption rate in the
	 * modified trap-mutation handler.
	 *
	 * @param i The position on the grid
	 * @return The rate
	 */
	double getLeftSideRate(int i) const override;

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
	std::vector<int> getConnectivity() const override;

	/**
	 * Tell reactant to output a representation of its reaction coefficients
	 * to the given output stream.
	 *
	 * @param os Output stream on which to output coefficients.
	 */
	virtual void outputCoefficientsTo(std::ostream& os) const override;

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
};

} /* end namespace xolotlCore */
#endif
