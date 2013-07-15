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
#include <PSIClusterReactionNetwork.h>
#include <string>

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
 * single-species He, V and I clusters first and all mixed clusters coming
 * last. Each species is ordered from the smallest cluster size, (1), to the
 * maximum size for that cluster. Instances of the appropriate cluster type are
 * instantiated during the loading process, but returned as PSIClusters.
 *
 * The ReactionNetwork's map of properties will contains the following
 * information about the network with the following keys:
 * > maxHeClusterSize - The number of He atoms in the largest single-species
 *  He cluster.
 * > maxVClusterSize - The number of atomic vacancies in the largest
 * single-species V cluster.
 * > maxIClusterSize - The number of interstitials in the largest
 * single-species I cluster.
 * > maxMixedClusterSize - The number of species of all types in the largest
 * mixed species in the network. It is equal to the sum of the max single
 * species helium and vacancy cluster sizes by default.
 * > numHeClusters - The number of single-species He clusters of all sizes in
 * the network.
 * > numVClusters - The number of single-species V clusters of all sizes in the
 * network.
 * > numIClusters - The number of single-species I clusters of all sizes in the
 * network.
 * > numMixedClusters - The number of mixed-species clusters of all sizes in the
 * network.
 */
class PSIClusterNetworkLoader {

private:

	/**
	 * The istream from which the network of clusters will be read.
	 */
	std::shared_ptr<std::istream> networkStream;

	/**
	 * The internal list of He clusters
	 */
	std::vector<std::shared_ptr<PSICluster>> heClusters;

	/**
	 * The internal list of vacancy clusters
	 */
	std::vector<std::shared_ptr<PSICluster>> vClusters;

	/**
	 * The internal list of interstitial clusters
	 */
	std::vector<std::shared_ptr<PSICluster>> iClusters;

	/**
	 * The internal list of HeV clusters
	 */
	std::vector<std::shared_ptr<PSICluster>> heVClusters;

	/**
	 * This operation creates a singles-species cluster of helium, vacancies or
	 * interstitials. It adds the cluster to the appropriate internal list of
	 * clusters for that type.
	 * @param numHe - The number of helium atoms
	 * @param numV - The number of atomic vacancies
	 * @param numI - The number of interstitial defects
	 * @return The new cluster
	 */
	std::shared_ptr<PSICluster> createCluster(int numHe, int numV, int numI,
			std::shared_ptr<std::map<std::string, std::string>> props);

public:

	/**
	 * The default constructor. The setInputstream() operation must be called
	 * if this constructor is used.
	 */
	PSIClusterNetworkLoader() {
	}

	/**
	 * An alternative constructor provided for convenience.
	 * @param inputstream The inputstream from which the cluster data should be
	 * loaded.
	 */
	PSIClusterNetworkLoader(const std::shared_ptr<std::istream> stream);

	/**
	 * Destructor
	 */
	virtual ~PSIClusterNetworkLoader() {
	}

	/**
	 * This operation specifies the inputstream from which cluster data should
	 * be loaded.
	 * @param inputstream The inputstream from which the cluster data should be
	 * loaded.
	 */
	void setInputstream(const std::shared_ptr<std::istream> stream);

	std::shared_ptr<std::istream> getInputstream();

	/**
	 * This operation will load the reaction network from the inputstream in
	 * the format specified previously. The network will be empty if it can not
	 * be loaded.
	 * @param network The reaction network
	 */
	std::shared_ptr<PSIClusterReactionNetwork> load();
};

} /* namespace xolotlCore */
#endif /* PSICLUSTERNETWORKLOADER_H_ */
