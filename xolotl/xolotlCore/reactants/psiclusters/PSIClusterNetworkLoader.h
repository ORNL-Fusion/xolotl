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

	// I formation energies in eV
	std::vector<double> iFormationEnergies = { 10.0, 18.5, 27.0, 35.0, 42.5,
			48.0 };
	// I diffusion factors in nm^2/s
	std::vector<double> iDiffusion = { 8.8e+10, 8.0e+10, 3.9e+10, 2.0e+10,
			1.0e+10 };
	// I migration energies in eV
	std::vector<double> iMigration = { 0.01, 0.02, 0.03, 0.04, 0.05 };

	// He formation energies in eV
	std::vector<double> heFormationEnergies = { 6.15, 11.44, 16.35, 21.0, 26.1,
			30.24, 34.93, 38.80 };
	// He diffusion factors in nm^2/s
	std::vector<double> heDiffusion = { 2.9e+10, 3.2e+10, 2.3e+10, 1.7e+10,
			5.0e+09, 1.0e+09, 5.0e+08 };
	// He migration energies in eV
	std::vector<double> heMigration = { 0.13, 0.20, 0.25, 0.20, 0.12, 0.3, 0.4 };

	// The diffusion factor for a single deuterium.
	double dOneDiffusionFactor = 2.83e+11;
	// The migration energy for a single deuterium.
	double dOneMigrationEnergy = 0.38;

	// The diffusion factor for a single tritium.
	double tOneDiffusionFactor = 2.31e+11;
	// The migration energy for a single tritium.
	double tOneMigrationEnergy = 0.38;

	// The diffusion factor for a single vacancy in nm^2/s
	double vOneDiffusion = 1.8e+12;
	// The migration energy for a single vacancy in eV
	double vOneMigration = 1.30;

	// The factor between He and H radius sizes
	double hydrogenRadiusFactor = 0.25;

	/**
	 * The maximum number of helium atoms that can be combined with a vacancy
	 * cluster with size equal to the index i in the array plus one. For
	 * example, an HeV size cluster with size 1 would have size = i+1 = 1 and i
	 * = 0. It could support a mixture of up to nine helium atoms with one
	 * vacancy.
	 */
	std::vector<int> maxHePerV = { 9, 14, 18, 20, 27, 30, 35, 40, 45, 50, 55,
			60, 65, 70, 75, 80, 85, 90, 95, 98, 100, 101, 103, 105, 107, 109,
			110, 112, 116 };

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
	 * This operation creates a super cluster from its list of cluster coordinates.
	 *
	 * @param list The list of coordinates composing this cluster
	 * @return The new cluster
	 */
	std::unique_ptr<PSICluster> createPSISuperCluster(std::set<std::tuple<int, int, int, int> > &list,
			IReactionNetwork& network) const;

	/**
	 * This operation creates a cluster.
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
