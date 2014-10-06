#ifndef VCLUSTER_H
#define VCLUSTER_H

// Includes
#include "PSICluster.h"
#include <xolotlPerf.h>

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
	VCluster() :
		PSICluster(1)
	{ }

public:

	/**
	 * The constructor. All VClusters must be initialized with a size.
	 * @param nV the number of atomic vacancies in the cluster
	 */
	VCluster(int nV, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	//! Destructor
	~VCluster();

	/**
	 * This operation returns a Reactant that is created using the copy
	 * constructor of VCluster.
	 * @return A copy of this reactant.
	 */
	virtual std::shared_ptr<Reactant> clone();

protected:
	/**
	 * Computes a row of the reaction connectivity matrix corresponding to
	 * this reactant.
	 *
	 * If two reactants alone can form a reaction, the element at the position
	 * of the second reactant is 1, otherwise 0.
	 */
	void createReactionConnectivity();

	/**
	 * Computes a row (or column) of the dissociation connectivity matrix
	 * corresponding to this cluster.
	 *
	 * Connections are made between this cluster and any clusters it affects
	 * in a dissociation reaction.
	 *
	 * The base-class implementation handles dissociation for regular clusters
	 * by processing the reaction.
	 *
	 */
	void createDissociationConnectivity();

};
//end class VCluster

} /* end namespace xolotlCore */

#endif
