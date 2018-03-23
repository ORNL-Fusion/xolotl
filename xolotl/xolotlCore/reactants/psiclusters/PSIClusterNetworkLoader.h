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
#include "PSIClusterReactionNetwork.h"

namespace xolotlCore {

/**
 * This class will load a reaction network composed of PSIClusters from an
 * inputstream.
 *
 * Lines of comments starting with a "#" will be ignored as will lines that do
 * not clearly provide the information above.
 */
class PSIClusterNetworkLoader: public NetworkLoader {

protected:

	/**
	 * The vacancy size at which the grouping scheme starts
	 */
	int vMin;

	/**
	 * The maximum size for helium clusters
	 */
	int maxHe;

	/**
	 * The maximum size for deuterium clusters
	 */
	int maxD;

	/**
	 * The maximum size for tritium clusters
	 */
	int maxT;

	/**
	 * The maximum size for interstitial clusters
	 */
	int maxI;

	/**
	 * The maximum size for vacancy clusters
	 */
	int maxV;

	/**
	 * The width of the group.
	 */
	int sectionWidth[4] = {};

	/**
	 * The list of clusters that will be grouped.
	 */
	std::set<std::tuple<int, int, int, int> > heVList;

	/**
	 * Private nullary constructor.
	 */
	PSIClusterNetworkLoader() :
			NetworkLoader(), vMin(1000000), maxHe(0), maxI(0), maxV(0), maxD(0), maxT(0) {
	}

	/**
	 * This operation creates a singles-species cluster of helium, vacancies or
	 * interstitials. It adds the cluster to the appropriate internal list of
	 * clusters for that type.
	 *
	 * @param numHe The number of helium atoms
	 * @param numD The number of deuterium atoms
	 * @param numT The number of tritium atoms
	 * @param numV The number of atomic vacancies
	 * @param numI The number of interstitial defects
	 * @return The new cluster
	 */
	std::unique_ptr<PSICluster> createPSICluster(int numHe, int numD, int numT, int numV, int numI,
			IReactionNetwork& network) const;

	/**
	 * This operation will add the given cluster to the network and reactants vector
	 * as a standard cluster or a dummy one if we do not want the reactions to happen.
	 *
	 * @param network The network
	 * @param reactants The vector of reactants kept by the loader
	 * @param cluster The cluster to add to them
	 */
	virtual void pushPSICluster(
			std::unique_ptr<PSIClusterReactionNetwork> & network,
			std::vector<std::reference_wrapper<Reactant> > & reactants,
			std::unique_ptr<PSICluster> & cluster);

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
	PSIClusterNetworkLoader(
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

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
	virtual ~PSIClusterNetworkLoader() {
	}

	/**
	 * This operation will load the reaction network from the inputstream in
	 * the format specified previously. The network will be empty if it can not
	 * be loaded.
	 *
	 * @return network The reaction network
	 */
	virtual std::unique_ptr<IReactionNetwork> load(const IOptions& options)
			override;

	/**
	 * This operation will generate the reaction network from options.
	 * The network will be empty if it can not be loaded.
	 *
	 * @param options The command line options
	 * @return network The reaction network
	 */
	virtual std::unique_ptr<IReactionNetwork> generate(const IOptions &options)
			override;

	/**
	 * This operation will apply a sectional grouping method to the network.
	 *
	 * @param The network to be modified.
	 */
	void applySectionalGrouping(PSIClusterReactionNetwork& network);

	/**
	 * This operation will set the helium size at which the grouping scheme starts.
	 *
	 * @param min The value for the size
	 */
	void setVMin(int min) {
		vMin = min;
	}

	/**
	 * This operation will set the helium width for the grouping scheme.
	 *
	 * @param w The value of the width
	 * @param axis The direction for the width
	 */
	void setWidth(int w, int axis) {
		sectionWidth[axis] = w;
	}
};

} /* namespace xolotlCore */

#endif /* PSICLUSTERNETWORKLOADER_H_ */
