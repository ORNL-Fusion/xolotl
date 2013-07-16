#ifndef INTERSTITIALCLUSTER_H
#define INTERSTITIALCLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of interstitial defects.
class InterstitialCluster: public PSICluster {

protected:

	/**
	 * This operation overrides getReactionConnectivity() from
	 * the base class to provide the proper reaction connectivity
	 * for an Interstitial cluster.
	 * @return The connectivity array
	 */
	std::vector<int> getReactionConnectivity();

	/**
	 * Computes a row of the dissociation connectivity matrix corresponding to
	 * this reactant.
	 *
	 * If two reactants together can be produced by a single reaction,
	 * the element at the position of the second reactant is 1, otherwise 0.
	 */
	virtual std::vector<int> getDissociationConnectivity();

public:

	/**
	 * The constructor. All InterstitialClusters must be initialized with
	 * a size.
	 */
	InterstitialCluster(int nI);

	//! Destructor
	~InterstitialCluster();

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
	 * Return whether or not this PSICluster is a product
	 * of the reaction between reactantI and reactantJ in
	 * this Reactants ReactionNetwork. This method should be
	 * specialized by subclasses to indicate whether or not they
	 * are the product of the given reaction.
	 *
	 * @param reactantI
	 * @param reactantJ
	 * @return
	 */
	virtual bool isProductReactant(int reactantI, int reactantJ);

	/**
	 * This virtual method is for subclasses to specialize
	 * and should return the reaction radius for the
	 * particular PSICluster subclass.
	 *
	 * @return
	 */
	virtual double getReactionRadius();

};
//end class InterstitialCluster

} /* end namespace xolotlCore */
#endif
