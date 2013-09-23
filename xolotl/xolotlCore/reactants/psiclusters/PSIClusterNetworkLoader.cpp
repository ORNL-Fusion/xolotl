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
#include <vector>
#include "HeCluster.h"
#include "VCluster.h"
#include "InterstitialCluster.h"
#include "HeVCluster.h"
// #include "HeInterstitialCluster.h"
#include "PSIClusterReactionNetwork.h"

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
		int numV, int numI) {

	// Local Declarations
	int clusterSize = 0;
	std::shared_ptr<PSICluster> cluster;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if (numHe > 0 && numV > 0) {
		// Create a new HeVCluster
		cluster = std::make_shared<HeVCluster>(numHe, numV);
		clusters.push_back(cluster);
	} else if (numHe > 0 && numI > 0) {
		throw std::string("HeliumInterstitialCluster is not yet implemented.");
		// FIXME! Add code to add it to the list
	} else if (numHe > 0) {
		// Create a new HeCluster
		cluster = std::make_shared<HeCluster>(numHe);
		clusters.push_back(cluster);
	} else if (numV > 0) {
		// Create a new VCluster
		cluster = std::make_shared<VCluster>(numV);
		clusters.push_back(cluster);
	} else if (numI > 0) {
		// Create a new ICluster
		cluster = std::make_shared<InterstitialCluster>(numI);
		// Add it to the ICluster list
		clusters.push_back(cluster);
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
std::shared_ptr<PSIClusterReactionNetwork> PSIClusterNetworkLoader::load() {

	// Local Declarations
	TokenizedLineReader<std::string> reader;
	std::vector<std::string> loadedLine;
	std::shared_ptr<PSIClusterReactionNetwork> network(
			new PSIClusterReactionNetwork());
	std::istringstream dataStream;
	std::string error(
			"PSIClusterNetworkLoader Exception: Insufficient or erroneous data.");
	int numHe = 0, numV = 0, numI = 0;
	double heBindingE = 0.0, vBindingE = 0.0, iBindingE = 0.0,
			trapMutationBindingE = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	bool mixed = false;
	std::vector<double> bindingEnergies;
	std::vector<std::shared_ptr<Reactant> > reactants;

	// Load the network if the stream is available
	if (networkStream != NULL) {
		// Load the stream
		reader.setInputStream(networkStream);

		// Loop over each line of the file, which should each be PSIClusters.
		loadedLine = reader.loadLine();
		while (loadedLine.size() > 0) {
			// Check the size of the loaded line
			if (loadedLine.size() < 9)
				// And notify the calling function if the size is insufficient
				throw error;
			// Load the sizes
			if (loadedLine[0][0] != '#') {
				numHe = std::stoi(loadedLine[0]);
				numV = std::stoi(loadedLine[1]);
				numI = std::stoi(loadedLine[2]);
				// Create the cluster
				auto nextCluster = createCluster(numHe, numV, numI);
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
				// Add the cluster to the network
				network->add(nextCluster);
				// Add it to the list so that we can set the network later
				reactants.push_back(nextCluster);
			}

			// Load the next line
			loadedLine = reader.loadLine();
		}

		for (auto reactantsIt = reactants.begin();
				reactantsIt != reactants.end(); ++reactantsIt) {
			(*reactantsIt)->setReactionNetwork(network);
		}

	}

	return network;
}
