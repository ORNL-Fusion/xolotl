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
#include "xolotlPerf/xolotlPerf.h"

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
	std::shared_ptr<PSICluster> cluster;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if (numHe > 0 && numV > 0) {
		// Create a new HeVCluster
		cluster = std::make_shared<HeVCluster>(numHe, numV, handlerRegistry);

		clusters.push_back(cluster);
	} else if (numHe > 0 && numI > 0) {
		throw std::string("HeliumInterstitialCluster is not yet implemented.");
		// FIXME! Add code to add it to the list
	} else if (numHe > 0) {
		// Create a new HeCluster
		cluster = std::make_shared<HeCluster>(numHe, handlerRegistry);

		clusters.push_back(cluster);
	} else if (numV > 0) {
		// Create a new VCluster
		cluster = std::make_shared<VCluster>(numV, handlerRegistry);

		clusters.push_back(cluster);
	} else if (numI > 0) {
		// Create a new ICluster
		cluster = std::make_shared<InterstitialCluster>(numI, handlerRegistry);

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
		const std::shared_ptr<std::istream> stream,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	handlerRegistry = registry;
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
	std::shared_ptr<PSIClusterReactionNetwork> network = std::make_shared<
			PSIClusterReactionNetwork>(handlerRegistry);

	std::string error(
			"PSIClusterNetworkLoader Exception: Insufficient or erroneous data.");
	int numHe = 0, numV = 0, numI = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::shared_ptr<Reactant> > reactants;

	// Load the network if the stream is available
	if (networkStream != NULL) {
		// Load the stream
		reader.setInputStream(networkStream);

		// Loop over each line of the file, which should each be PSIClusters.
		loadedLine = reader.loadLine();
		while (loadedLine.size() > 0) {
			// Check the size of the loaded line
			if (loadedLine.size() < 6)
				// And notify the calling function if the size is insufficient
				throw error;
			// Load the sizes
			if (loadedLine[0][0] != '#') {
				numHe = std::stoi(loadedLine[0]);
				numV = std::stoi(loadedLine[1]);
				numI = std::stoi(loadedLine[2]);
				// Create the cluster
				auto nextCluster = createCluster(numHe, numV, numI);
				// Load the energies
				formationEnergy = convertStrToDouble(loadedLine[3]);
				migrationEnergy = convertStrToDouble(loadedLine[4]);
				diffusionFactor = convertStrToDouble(loadedLine[5]);
				// Set the formation energy
				nextCluster->setFormationEnergy(formationEnergy);
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

		// Set the network for all of the reactants. This MUST be done manually.
		for (auto reactantsIt = reactants.begin();
				reactantsIt != reactants.end(); ++reactantsIt) {
			(*reactantsIt)->setReactionNetwork(network);
		}

	}

	return network;
}
