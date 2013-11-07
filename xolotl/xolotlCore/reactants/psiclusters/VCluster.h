#ifndef VCLUSTER_H
#define VCLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of atomic vacancies.
 */
class VCluster: public PSICluster {

private:

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	VCluster() : PSICluster(1) {}

public:

	/**
	 * The constructor. All VClusters must be initialized with a size.
	 * @param nV the number of atomic vacancies in the cluster
	 */
	VCluster(int nV);

	//! Destructor
	~VCluster();

	/**
	 * This operation returns a Reactant that is created using the copy
	 * constructor of VCluster.
	 * @return A copy of this reactant.
	 */
	virtual std::shared_ptr<Reactant> clone();

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
	virtual bool isProductReactant(const Reactant & reactantI, const Reactant & reactantJ);

	/**
	 * This virtual method is for subclasses to specialize
	 * and should return the reaction radius for the
	 * particular PSICluster subclass.
	 *
	 * @return
	 */
	virtual double getReactionRadius() const;
	
	/**
	 * This operation returns the compositon of this reactant.
	 * @return The composition returned as a map with keys naming distinct
	 * elements and values indicating the amount of the element present.
	 */
	virtual std::map<std::string,int> getComposition() const;

protected:
	/**
	 * Computes a row of the reaction connectivity matrix corresponding to
	 * this reactant.
	 *
	 * If two reactants alone can form a reaction, the element at the position
	 * of the second reactant is 1, otherwise 0.
	 */
	void createReactionConnectivity();

};
//end class VCluster

} /* end namespace xolotlCore */

#endif
