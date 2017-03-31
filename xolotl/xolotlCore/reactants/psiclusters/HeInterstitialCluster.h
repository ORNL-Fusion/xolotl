#ifndef HEINTERSTITIALCLUSTER_H
#define HEINTERSTITIALCLUSTER_H

// Includes
#include "PSICluster.h"
#include <string>
#include <map>

namespace xolotlCore {

/**
 *  A cluster composed of helium and intersititial.
 */
class HeInterstitialCluster : public PSICluster {

private:

	//! The number of helium atoms in this cluster.
	int numHe;

	//! The number of interstitial defects in this cluster.
	int numI;

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	HeInterstitialCluster() :
		PSICluster()
	{ numHe = 0; numI = 0; }

public:

	/**
	 * The constructor. All HeInterstitialClusters must be initialized with a map
	 * that describes the species of which the cluster is composed. The map
	 * should contain as its keys the names of the species and the sizes of the
	 * species as its values. The names of the species must be one of
	 * {He,I}.
	 *
	 * @param numHe The number of helium atoms in this cluster
	 * @param numI The number of interstitial defect in this cluster
	 * @param registry The performance handler registry
	 */
	HeInterstitialCluster(int numHe, int numI,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor
	 *
	 * @param other the reactant to be copied
	 */
	HeInterstitialCluster(HeInterstitialCluster &other);

	//! Destructor
	~HeInterstitialCluster() {}

	/**
	 * Returns a reactant created using the copy constructor
	 */
	virtual std::shared_ptr<IReactant> clone() {
		return std::shared_ptr<IReactant> (new HeInterstitialCluster(*this));
	}

	/**
	 * This operation returns true to signify that this cluster is a mixture of
	 * He and I.
	 *
	 * @return True if mixed
	 */
	virtual bool isMixed() const {return true;}
};
//end class HeInterstitialCluster

} /* end namespace xolotlCore */
#endif
