#ifndef NECLUSTERNETWORKLOADER_H_
#define NECLUSTERNETWORKLOADER_H_

//Includes
#include <NECluster.h>
#include <NetworkLoader.h>

namespace xolotlCore {

/**
 * This class will load a reaction network composed of NEClusters from an
 * HDF5 file.
 *
 * The network will be returned as a ReactionNetwork of NEClusters ordered with
 * single-species Xe, V and I clusters first and all mixed clusters coming
 * last. Each species is ordered from the smallest cluster size, (1), to the
 * maximum size for that cluster. Instances of the appropriate cluster type are
 * instantiated during the loading process, but returned as NEClusters.
 *
 * The ReactionNetwork's map of properties will contains the following
 * information about the network with the following keys:
 * > maxXeClusterSize - The number of Xe atoms in the largest single-species
 *  Xe cluster.
 * > maxVClusterSize - The number of atomic vacancies in the largest
 * single-species V cluster.
 * > maxIClusterSize - The number of interstitials in the largest
 * single-species I cluster.
 * > maxMixedClusterSize - The number of species of all types in the largest
 * mixed species in the network. It is equal to the sum of the max single
 * species helium and vacancy cluster sizes by default.
 * > numXeClusters - The number of single-species Xe clusters of all sizes in
 * the network.
 * > numVClusters - The number of single-species V clusters of all sizes in the
 * network.
 * > numIClusters - The number of single-species I clusters of all sizes in the
 * network.
 * > numMixedClusters - The number of mixed-species clusters of all sizes in the
 * network.
 */
class NEClusterNetworkLoader : public NetworkLoader {

protected:

	/**
	 * The xenon size at which the grouping scheme starts
	 */
	int xeMin;

	/**
	 * The width of the group in the xenon direction.
	 */
	int sectionWidth;

	/**
	 * Private nullary constructor.
	 */
	NEClusterNetworkLoader() {}

	/**
	 * This operation creates a cluster of helium, vacancies and/or
	 * interstitials. It adds the cluster to the appropriate internal list of
	 * clusters for that type.
	 *
	 * @param numXe The number of helium atoms
	 * @param numV The number of atomic vacancies
	 * @param numI The number of interstitial defects
	 * @return The new cluster
	 */
	std::shared_ptr<NECluster> createNECluster(int numXe, int numV, int numI);

public:

	/**
	 * The default constructor. The setInputstream() operation must be called
	 * if this constructor is used.
	 *
	 * @param registry The performance handler registry
	 */
	NEClusterNetworkLoader(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * An alternative constructor provided for convenience.
	 *
	 * @param stream The inputstream from which the cluster data should be
	 * loaded
	 * @param registry The performance handler registry
	 */
	NEClusterNetworkLoader(const std::shared_ptr<std::istream> stream,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Destructor
	 */
	virtual ~NEClusterNetworkLoader() {}

	/**
	 * This operation will load the reaction network from the inputstream in
	 * the format specified previously. The network will be empty if it can not
	 * be loaded.
	 *
	 * @return network The reaction network
	 */
	virtual std::shared_ptr<IReactionNetwork> load();

	/**
	 * This operation will generate the reaction network from options.
	 * The network will be empty if it can not be loaded.
	 *
	 * @param options The command line options
	 * @return network The reaction network
	 */
	virtual std::shared_ptr<IReactionNetwork> generate(IOptions &options);

	/**
	 * This operation will apply a grouping method to the network.
	 *
	 * @param The network to be modified.
	 */
	void applyGrouping(std::shared_ptr<IReactionNetwork> network);

	/**
	 * This operation will set the xenon size at which the grouping scheme starts.
	 *
	 * @param min The value for the size
	 */
	void setXeMin (int min) {xeMin = min;}

	/**
	 * This operation will set the xenon width for the grouping scheme.
	 *
	 * @param w The value of the width
	 */
	void setWidth (int w) {sectionWidth = w;}
};

} /* namespace xolotlCore */

#endif /* NECLUSTERNETWORKLOADER_H_ */
