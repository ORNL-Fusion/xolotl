/*
 * PSIClusterNetworkLoader.h
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#ifndef PSICLUSTERNETWORKLOADER_H_
#define PSICLUSTERNETWORKLOADER_H_

//Includes
#include <PSICluster.h>
#include <NetworkLoader.h>

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
 * > The formation energy
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
class PSIClusterNetworkLoader : public NetworkLoader {

protected:

	/**
	 * The vacancy size at which the grouping scheme starts
	 */
	int vMin;

	/**
	 * The width of the group in the helium direction.
	 */
	int heSectionWidth;

	/**
	 * The width of the group in the vacancy direction.
	 */
	int vSectionWidth;

	/**
	 * Private nullary constructor.
	 */
	PSIClusterNetworkLoader() {}

	/**
	 * This operation creates a singles-species cluster of helium, vacancies or
	 * interstitials. It adds the cluster to the appropriate internal list of
	 * clusters for that type.
	 *
	 * @param numHe The number of helium atoms
	 * @param numV The number of atomic vacancies
	 * @param numI The number of interstitial defects
	 * @return The new cluster
	 */
	std::shared_ptr<PSICluster> createPSICluster(int numHe, int numV, int numI);

	/**
	 * This operation computes the formation energy associated to the
	 * cluster of the given size.
	 *
	 * @param numHe The number of helium atoms
	 * @param numV The number of atomic vacancies
	 * @return The corresponding formation energy
	 */
	double getHeVFormationEnergy(int numHe, int numV);

public:

	/**
	 * The default constructor. The setInputstream() operation must be called
	 * if this constructor is used.
	 *
	 * @param registry The performance handler registry
	 */
	PSIClusterNetworkLoader(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * An alternative constructor provided for convenience.
	 *
	 * @param stream The inputstream from which the cluster data should be
	 * loaded
	 * @param registry The performance handler registry
	 */
	PSIClusterNetworkLoader(const std::shared_ptr<std::istream> stream,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Destructor
	 */
	virtual ~PSIClusterNetworkLoader() {}

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
	 * This operation will apply a sectional grouping method to the network.
	 *
	 * @param The network to be modified.
	 */
	void applySectionalGrouping(std::shared_ptr<IReactionNetwork> network);

	/**
	 * This operation will set the helium size at which the grouping scheme starts.
	 *
	 * @param min The value for the size
	 */
	void setVMin (int min) {vMin = min;}

	/**
	 * This operation will set the helium width for the grouping scheme.
	 *
	 * @param w The value of the width
	 */
	void setHeWidth (int w) {heSectionWidth = w;}

	/**
	 * This operation will set the vacancy width for the grouping scheme.
	 *
	 * @param w The value of the width
	 */
	void setVWidth (int w) {vSectionWidth = w;}
};

} /* namespace xolotlCore */

#endif /* PSICLUSTERNETWORKLOADER_H_ */
