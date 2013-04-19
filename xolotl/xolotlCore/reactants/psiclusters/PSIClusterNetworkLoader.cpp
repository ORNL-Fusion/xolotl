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
#include "HeCluster.h"

using namespace xolotlCore;

/**
 * This operation converts a string to a double, taking in to account the fact
 * that the input file may contain keys such as "infinite."
 * @param inString the string to be converted
 * @return the string as a double
 */
static inline double convertStrToDouble(const std::string inString) {
	// Conv
	return (inString.compare("infinite") == 0) ?
			std::numeric_limits<double>::infinity() :
			strtod(inString.c_str(), NULL);
}

/**
 * This operation creates a singles-species cluster of helium, vacancies or interstitials.
 * @param numHe - The number of helium atoms
 * @param numV - The number of atomic vacancies
 * @param numI - The number of interstitial defects
 * @return The new cluster
 */
static PSICluster createSingleSpeciesCluster(int numHe, int numV, int numI) {
	// Local Declarations

	HeCluster cluster(1);

	return cluster;
}

/**
 * This operation creates a singles-species cluster of helium, vacancies or interstitials.
 * @param numHe - The number of helium atoms
 * @param numV - The number of atomic vacancies
 * @param numI - The number of interstitial defects
 * @return The new cluster
 */
static PSICluster createMixedSpeciesCluster(int numHe, int numV, int numI) {
	PSICluster cluster(1);

	return cluster;
}

/**
 * An alternative constructor provided for convenience.
 * @param inputstream The inputstream from which the cluster data should be
 * loaded.
 */
PSIClusterNetworkLoader::PSIClusterNetworkLoader(
		const std::shared_ptr<std::istream> stream) {
}

/**
 * This operation specifies the inputstream from which cluster data should
 * be loaded.
 * @param inputstream The inputstream from which the cluster data should be
 * loaded.
 */
void PSIClusterNetworkLoader::setInputstream(
		const std::shared_ptr<std::istream> stream) {

	// Keep the stream, but only if it is not NULL
	if (stream)
		networkStream = stream;

	return;
}

/**
 * This operation will load the reaction network from the inputstream in
 * the format specified previously. The network will be empty if it can not
 * be loaded.
 * @param network The reaction network
 */
std::shared_ptr<ReactionNetwork> PSIClusterNetworkLoader::load() {

	// Local Declarations
	TokenizedLineReader<std::string> reader;
	std::vector<std::string> loadedLine;
	std::shared_ptr<ReactionNetwork> network(new ReactionNetwork());
	std::istringstream dataStream;
	int numHe = 0, numV = 0, numI = 0;
	double heBindingE = 0.0, vBindingE = 0.0, iBindingE = 0.0,
			trapMutationBindingE = 0.0;
	bool mixed = false;

	// Load the network if the stream is available
	if (networkStream != NULL) {
		// Load the stream
		reader.setInputStream(networkStream);
		// Loop over each line of the file, which should each be PSIClusters.
		loadedLine = reader.loadLine();
		while (loadedLine.size() > 0) {
			// Load the sizes
			numHe = strtol(loadedLine.at(0).c_str(), NULL, 10);
			numV = strtol(loadedLine.at(1).c_str(), NULL, 10);
			numI = strtol(loadedLine.at(2).c_str(), NULL, 10);
			// Load the binding energies
			heBindingE = convertStrToDouble(loadedLine.at(3));
			vBindingE = convertStrToDouble(loadedLine.at(4));
			iBindingE = convertStrToDouble(loadedLine.at(5));
			trapMutationBindingE = convertStrToDouble(loadedLine.at(6));
			// Determine whether or not this is a mixed cluster
			mixed = (((numHe > 0) + (numV > 0) + (numI > 0)) <= 1);
			// Create the cluster
			PSICluster nextCluster =
					(!mixed) ?
							createSingleSpeciesCluster(numHe,numV,numI) :
							createMixedSpeciesCluster(numHe,numV,numI);
			// Load the next line
			loadedLine = reader.loadLine();
		}
	}

	return network;

}
