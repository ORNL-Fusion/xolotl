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


}; //end class XeCluster

} /* namespace xolotlCore */
#endif
