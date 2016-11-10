#ifndef XECLUSTER_H
#define XECLUSTER_H

// Includes
#include "NECluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of helium.
 */
class XeCluster: public NECluster {

private:

	/**
	 * The default constructor is private because NEClusters must always be
	 * initialized with a size and performance handler registry
	 */
	XeCluster() :
		NECluster() {}

public:

	/**
	 * The constructor. All XeClusters must be initialized with a size.
	 *
	 * @param nXe the number of helium atoms in the cluster
	 * @param registry The performance handler registry
	 */
	XeCluster(int nXe, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Destructor
	 */
	~XeCluster() {}

	/**
	 * Returns a reactant created using the copy constructor
	 */
	virtual std::shared_ptr<IReactant> clone() {
		return std::shared_ptr<IReactant> (new XeCluster(*this));
	}

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
	 * Computes a row of the dissociation connectivity matrix
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


}; //end class XeCluster

} /* namespace xolotlCore */
#endif
