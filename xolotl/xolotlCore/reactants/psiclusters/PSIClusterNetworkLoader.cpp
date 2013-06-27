/*
 * PSIClusterNetworkLoader.cpp
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#include "PSIClusterNetworkLoader.h"
#include <TokenizedLineReader.h>
#include <stdio.h>
#include <limits>
#include <algorithm>
#include "HeCluster.h"
#include "VCluster.h"
#include "InterstitialCluster.h"
#include "MixedSpeciesCluster.h"

using namespace xolotlCore;

/**
 * This operation converts a string to a double, taking in to account the fact
 * that the input file may contain keys such as "infinite."
 * @param inString the string to be converted
 * @return the string as a double
 */
static inline double convertStrToDouble(const std::string inString) {
	return (inString.compare("infinite") == 0) ?
			std::numeric_limits<double>::infinity() :
			strtod(inString.c_str(), NULL);
}

/**
 * This operation creates a singles-species cluster of helium, vacancies or
 * interstitials. It adds the cluster to the appropriate internal list of
 * clusters for that type.
 * @param numHe - The number of helium atoms
 * @param numV - The number of atomic vacancies
 * @param numI - The number of interstitial defects
 * @return The new cluster
 */
std::shared_ptr<PSICluster> PSIClusterNetworkLoader::createCluster(int numHe,
		int numV, int numI,
		std::shared_ptr<std::map<std::string, std::string>> props) {

	// Local Declarations
	std::map<std::string, int> speciesMap;
	int numClusters = 0;
	long clusterSize = 0, maxClusterSize = 0;
	std::string numClustersTag, maxClustersTag;
	std::shared_ptr<PSICluster> cluster;

	// Determine whether or not this is a mixed cluster
	bool mixed = (((numHe > 0) + (numV > 0) + (numI > 0)) > 1);

	/* Most of the clusters will be mixed - there are only about 30
	 * single-species clusters.
	 */
	if (mixed) {
		// Load the species map so that the mixed cluster can be properly created.
		speciesMap["He"] = numHe;
		speciesMap["V"] = numV;
		speciesMap["I"] = numI;
		cluster = std::make_shared < MixedSpeciesCluster > (speciesMap);
		numClustersTag = "numMixedClusters";
		// Update the number of clusters
		numClusters = strtol((*props)[numClustersTag].c_str(), NULL, 10) + 1;
		(*props)[numClustersTag] = std::to_string(
				static_cast<long long>(numClusters));
		mixedClusters.push_back(cluster);
	} else {
		/* Switch over the three types, create the cluster and set the properties.
		 * Start with He as they are probably listed first.
		 */
		if (numHe > 0) {
			cluster = std::make_shared < HeCluster > (numHe);
			clusterSize = numHe;
			numClustersTag = "numHeClusters";
			maxClustersTag = "maxHeClusterSize";
			heClusters.push_back(cluster);
		} else if (numV > 0) { // Vacancies
			cluster = std::make_shared < VCluster > (numV);
			clusterSize = numV;
			numClustersTag = "numVClusters";
			maxClustersTag = "maxVClusterSize";
			vClusters.push_back(cluster);
		} else { // Default to interstitial defects.
			cluster = std::make_shared < InterstitialCluster > (numI);
			clusterSize = numI;
			numClustersTag = "numIClusters";
			maxClustersTag = "maxIClusterSize";
			iClusters.push_back(cluster);
		}
		// Update the number of clusters
		numClusters = strtol((*props)[numClustersTag].c_str(), NULL, 10) + 1;
		(*props)[numClustersTag] = std::to_string(
				static_cast<long long>(numClusters));
		// Update the max size if required - compute the max from the old and current values
		maxClusterSize = strtol((*props)[maxClustersTag].c_str(), NULL, 10);
		maxClusterSize = std::max(clusterSize, maxClusterSize);
		(*props)[maxClustersTag] = std::to_string(
				static_cast<long long>(clusterSize));
	}

	return cluster;
}

/**
 * An alternative constructor provided for convenience.
 * @param inputstream The inputstream from which the cluster data should be
 * loaded.
 */
PSIClusterNetworkLoader::PSIClusterNetworkLoader(
		const std::shared_ptr<std::istream> stream) {
	setInputstream(stream);
}

/**
 * This operation specifies the inputstream from which cluster data should
 * be loaded.
 * @param inputstream The inputstream from which the cluster data should be
 * loaded.
 */
void PSIClusterNetworkLoader::setInputstream(
		const std::shared_ptr<std::istream> stream) {

	networkStream = stream;
}

std::shared_ptr<std::istream> PSIClusterNetworkLoader::getInputstream() {
	return networkStream;
}

/**
 * This operation will load the reaction network from the inputstream in
 * the format specified previously. The network will be empty if it can not
 * be loaded.
 * @param network The reaction network
 *
 * This operation throws a std::string as an exception if there is a problem.
 */
std::shared_ptr<ReactionNetwork> PSIClusterNetworkLoader::load() {

	// Local Declarations
	TokenizedLineReader<std::string> reader;
	std::vector<std::string> loadedLine;
	std::shared_ptr<ReactionNetwork> network(new ReactionNetwork());
	std::istringstream dataStream;
	std::string error(
			"PSIClusterNetworkLoader Exception: Insufficient or erroneous data.");
	std::shared_ptr<std::map<std::string, std::string>> props;
	int numHe = 0, numV = 0, numI = 0;
	double heBindingE = 0.0, vBindingE = 0.0, iBindingE = 0.0,
			trapMutationBindingE = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	bool mixed = false;
	std::vector<double> bindingEnergies;

	// Load the network if the stream is available
	if (networkStream != NULL) {
		// Load the stream
		reader.setInputStream(networkStream);
		// Loop over each line of the file, which should each be PSIClusters.
		loadedLine = reader.loadLine();
		// Setup the properties map
		props = network->properties;
		(*props)["maxHeClusterSize"] = "0";
		(*props)["maxVClusterSize"] = "0";
		(*props)["maxIClusterSize"] = "0";
		(*props)["numHeClusters"] = "0";
		(*props)["numVClusters"] = "0";
		(*props)["numIClusters"] = "0";
		(*props)["numMixedClusters"] = "0";
		while (loadedLine.size() > 0) {
			// Check the size of the loaded line
			if (loadedLine.size() < 9)
				// And notify the calling function if the size is insufficient
				throw error;
			// Load the sizes
			numHe = strtol(loadedLine[0].c_str(), NULL, 10);
			numV = strtol(loadedLine[1].c_str(), NULL, 10);
			numI = strtol(loadedLine[2].c_str(), NULL, 10);
			// Create the cluster
			std::shared_ptr<PSICluster> nextCluster = createCluster(numHe, numV,
					numI, props);
			// Load the binding energies
			heBindingE = convertStrToDouble(loadedLine[3]);
			vBindingE = convertStrToDouble(loadedLine[4]);
			iBindingE = convertStrToDouble(loadedLine[5]);
			trapMutationBindingE = convertStrToDouble(loadedLine[6]);
			migrationEnergy = convertStrToDouble(loadedLine[7]);
			diffusionFactor = convertStrToDouble(loadedLine[8]);
			// Create the binding energies array and set it
			bindingEnergies.clear();
			bindingEnergies.push_back(heBindingE);
			bindingEnergies.push_back(vBindingE);
			bindingEnergies.push_back(iBindingE);
			bindingEnergies.push_back(trapMutationBindingE);
			nextCluster->setBindingEnergies(bindingEnergies);
			// Set the diffusion factor and migration energy
			nextCluster->setMigrationEnergy(migrationEnergy);
			nextCluster->setDiffusionFactor(diffusionFactor);
			// Load the next line
			loadedLine = reader.loadLine();
		}
		// Load the clusters into the network, starting with He
		for (int i = 0; i < heClusters.size(); i++)
			network->reactants->push_back(heClusters[i]);
		// Load the vacancies into the network
		for (int i = 0; i < vClusters.size(); i++)
			network->reactants->push_back(vClusters[i]);
		// Load the interstitials into the network
		for (int i = 0; i < iClusters.size(); i++)
			network->reactants->push_back(iClusters[i]);
		// Load the mixed species clusters into the network
		for (int i = 0; i < mixedClusters.size(); i++)
			network->reactants->push_back(mixedClusters[i]);
		// Set the max number of mixed clusters
		int maxMixedClusterSize = strtol(
				(*props)["maxHeClusterSize"].c_str(), NULL, 10)
				+ strtol((*props)["maxVClusterSize"].c_str(),
						NULL, 10)
				+ strtol((*props)["maxIClusterSize"].c_str(),
						NULL, 10);
		(*props)["maxMixedClusterSize"] = std::to_string((long long) maxMixedClusterSize);
		std::cout << "MSize = " << (*props)["maxMixedClusterSize"] << std::endl;

	}

	return network;

}
