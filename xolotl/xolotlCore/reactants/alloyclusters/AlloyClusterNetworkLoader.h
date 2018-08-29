#ifndef ALLOYCLUSTERNETWORKLOADER_H_
#define ALLOYCLUSTERNETWORKLOADER_H_

//Includes
#include <AlloyCluster.h>
#include <NetworkLoader.h>
#include <AlloyClusterReactionNetwork.h>

namespace xolotlCore {

/**
 * This class will load a reaction network composed of AlloyClusters from an
 * HDF5 file.
 *
 * The network will be returned as a ReactionNetwork of AlloyClusters.
 * Each species is ordered from the smallest cluster size, (1), to the
 * maximum size for that cluster. Instances of the appropriate cluster type are
 * instantiated during the loading process, but returned as AlloyClusters.
 */
class AlloyClusterNetworkLoader: public NetworkLoader {

protected:

	/**
	 * The size at which the grouping scheme starts
	 */
	int sizeMin;

	/**
	 * The width of the group.
	 */
	int sectionWidth;

	/**
	 * Private nullary constructor.
	 */
	AlloyClusterNetworkLoader() :
			sizeMin(-1), sectionWidth(-1) {
	}

	/**
	 * This operation creates a singles-species cluster of vacancies,
	 * interstitials, void, faulted, frank, perfect.
	 * It adds the cluster to the appropriate internal list of
	 * clusters for that type.
	 *
	 * @param numV The number of atomic vacancies
	 * @param numI The number of interstitial defects
	 * @param numVoid The number of void
	 * @param numFaulted The number of faulted
	 * @param numFrank The number of frank
	 * @param numPerfect The number of perfect
	 * @return The new cluster
	 */
	std::unique_ptr<AlloyCluster> createAlloyCluster(int numV, int numI,
			int numVoid, int numFaulted, int numFrank, int numPerfect,
			IReactionNetwork& network);

	/**
	 * This operation will add the given cluster to the network and reactants vector
	 * as a standard cluster or a dummy one if we do not want the reactions to happen.
	 *
	 * @param network The network
	 * @param reactants The vector of reactants kept by the loader
	 * @param cluster The cluster to add to them
	 */
	virtual void pushAlloyCluster(
			std::unique_ptr<AlloyClusterReactionNetwork> & network,
			std::vector<std::reference_wrapper<Reactant> > & reactants,
			std::unique_ptr<AlloyCluster> & cluster);

public:

	/**
	 * The default constructor. The setInputstream() operation must be called
	 * if this constructor is used.
	 *
	 * @param registry The performance handler registry
	 */
	AlloyClusterNetworkLoader(
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * An alternative constructor provided for convenience.
	 *
	 * @param stream The inputstream from which the cluster data should be
	 * loaded
	 * @param registry The performance handler registry
	 */
	AlloyClusterNetworkLoader(const std::shared_ptr<std::istream> stream,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Destructor
	 */
	virtual ~AlloyClusterNetworkLoader() {
	}

	/**
	 * This operation will load the reaction network from the inputstream in
	 * the format specified previously. The network will be empty if it can not
	 * be loaded.
	 *
	 * @param options The command line options
	 * @return network The reaction network
	 */
	virtual std::unique_ptr<IReactionNetwork> load(const IOptions &options);

	/**
	 * This operation will generate the reaction network from options.
	 * The network will be empty if it can not be loaded.
	 *
	 * @param options The command line options
	 * @return network The reaction network
	 */
	virtual std::unique_ptr<IReactionNetwork> generate(const IOptions &options);

	/**
	 * This operation will apply a grouping method to the network.
	 *
	 * @param The network to be modified.
	 */
	void applyGrouping(IReactionNetwork& network);

	/**
	 * This operation will set the size at which the grouping scheme starts.
	 *
	 * @param min The value for the size
	 */
	void setMin(int min) {
		sizeMin = min;
	}

	/**
	 * This operation will set the width for the grouping scheme.
	 *
	 * @param w The value of the width
	 */
	void setWidth(int w) {
		sectionWidth = w;
	}
};

} /* namespace xolotlCore */

#endif /* ALLOYCLUSTERNETWORKLOADER_H_ */
