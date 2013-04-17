/*
 * PSIClusterNetworkLoader.h
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#ifndef PSICLUSTERNETWORKLOADER_H_
#define PSICLUSTERNETWORKLOADER_H_

//Includes
#include <map>
#include <memory>
#include "PSICluster.h"
#include <ReactionNetwork.h>

namespace xolotlCore {

/**
 * This class will load a reaction network composed of PSIClusters from an
 * inputstream.
 *
 * The data in the stream should contain the information for a single cluster on
 * each line with the following quantities specified and separated by a single
 * space each:
 * > The number of He in the cluster
 * > The number of V in the cluster
 * > The number of I in the cluster
 * > The He binding energy
 * > The V binding energy
 * > The I binding energy
 * > The trap mutation binding energy
 *
 * Where appropriate, any binding energy may be specified as "Infinity" to
 * signify that the cluster does not undergo that type of dissociation.
 *
 * Lines of comments starting with a "#" will be ignored as will lines that do
 * not clearly provide the information above.
 *
 * The network will be returned as a ReactionNetwork of PSIClusters ordered with
 * single-species He, I and V clusters first and all mixed clusters coming
 * last. Instances of the appropriate cluster type are instantiated during the
 * loading process, but returned as instances of the base class.
 *
 * The ReactionNetwork's map of properties will contains the following
 * information about the network with the following keys:
 * > maxHeClusterSize - The number of He atoms in the largest single-species
 *  He cluster.
 * > maxVClusterSize - The number of atomic vacancies in the largest
 * single-species V cluster.
 * > maxInterstitialClusterSize - The number of interstitials in the largest
 * single-species I cluster.
 * > numHeClusters - The number of single-species He clusters of all sizes in
 * the network.
 * > numVClusters - The number of single-species V clusters of all sizes in the
 * network.
 * > numInterstitialClusters - The number of single-species I clusters of all
 * sizes in the network.
 */
class PSIClusterNetworkLoader {

private:

	/**
	 * The istream from which the network of clusters will be read.
	 */
	std::shared_ptr<std::istream> networkStream;

public:

	/**
	 * The default constructor. The setInputstream() operation must be called
	 * if this constructor is used.
	 */
	PSIClusterNetworkLoader();

	/**
	 * An alternative constructor provided for convenience.
	 * @param inputstream The inputstream from which the cluster data should be
	 * loaded.
	 */
	PSIClusterNetworkLoader(std::shared_ptr<std::istream> stream);

	/**
	 * Destructor
	 */
	virtual ~PSIClusterNetworkLoader();

	/**
	 * This operation specifies the inputstream from which cluster data should
	 * be loaded.
	 * @param inputstream The inputstream from which the cluster data should be
	 * loaded.
	 */
	void setInputstream(std::shared_ptr<std::istream> stream);

	/**
	 * This operation will load the reaction network from the inputstream in
	 * the format specified previously. The network will be empty if it can not
	 * be loaded.
	 * @param network The reaction network
	 */
	std::shared_ptr<ReactionNetwork> load();

};

} /* namespace xolotlCore */
#endif /* PSICLUSTERNETWORKLOADER_H_ */
