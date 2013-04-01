/*
 * PSIClusterNetworkLoader.h
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#ifndef PSICLUSTERNETWORKLOADER_H_
#define PSICLUSTERNETWORKLOADER_H_

//Includes
#include <stdlib.h>
#include <psicluster.h>

namespace xolotlCore {

/**
 * This class will load a reaction network composed of PSIClusters from an
 * inputstream.
 *
 * The data in the stream should contain the information for a single cluster on
 * each line with the following quantities specified and separated by a single
 * space each:
 * > The number of He in the cluster
 * > The number of I in the cluster
 * > The number of V in the cluster
 * > The He binding energy
 * > The I binding energy
 * > The V binding energy
 * > The trap mutation binding energy
 *
 * Lines of comments starting with a "#" will be ignored as will lines that do
 * not clearly provide the information above.
 *
 * The network will be returned as a std::vector of PSIClusters ordered with
 * single-species He, I and V clusters first and all mixed clusters coming
 * last. Instances of the appropriate cluster type are instantiated during the
 * loading process, but returned as instances of the base class.
 *
 * A map of properties will be produced that contains information about the
 * network with the following keys:
 * > maxHeClusterSize - The number of He atoms in the largest single-species
 *  He cluster.
 * > maxInterstitialClusterSize - The number of interstitials in the largest
 * single-species I cluster.
 * > maxIVClusterSize - The number of atomic vacancies in the largest
 * single-species V cluster.
 * > numHeClusters - The number of single-species He clusters of all sizes in
 * the network.
 * > numInterstitialClusters - The number of single-species I clusters of all
 * sizes in the network.
 * > numVClusters - The number of single-species V clusters of all sizes in the
 * network.
 *
 */
class PSIClusterNetworkLoader {

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
	PSIClusterNetworkLoader(std::istream inputstream);

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
	void setInputstream(std::istream inputstream);

	/**
	 * This operation will load the reaction network from the inputstream in
	 * the format specified previously. The network will be empty if it can not
	 * be loaded.
	 * @return The reaction network //Should it be vector<vector>?
	 */
	std::vector<PSICluster> load();

	/**
	 * This operation returns the properties map that was created for the
	 * network. If it is called before the network is loaded or the network can
	 * not be loaded, it will be empty.
	 * @return The property map
	 */
	std::map<std::string,std::string> getProperties();

};

} /* namespace xolotlCore */
#endif /* PSICLUSTERNETWORKLOADER_H_ */
