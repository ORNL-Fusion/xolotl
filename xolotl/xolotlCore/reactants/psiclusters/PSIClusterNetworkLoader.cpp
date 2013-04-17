/*
 * PSIClusterNetworkLoader.cpp
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#include "PSIClusterNetworkLoader.h"
#include <TokenizedLineReader.h>

using namespace xolotlCore;

/**
 * The default constructor. The setInputstream() operation must be called
 * if this constructor is used.
 */
PSIClusterNetworkLoader::PSIClusterNetworkLoader() {
}

/**
 * An alternative constructor provided for convenience.
 * @param inputstream The inputstream from which the cluster data should be
 * loaded.
 */
PSIClusterNetworkLoader::PSIClusterNetworkLoader(std::shared_ptr<std::istream> stream) {
}

/**
 * Destructor
 */
PSIClusterNetworkLoader::~PSIClusterNetworkLoader() {
}

/**
 * This operation specifies the inputstream from which cluster data should
 * be loaded.
 * @param inputstream The inputstream from which the cluster data should be
 * loaded.
 */
void PSIClusterNetworkLoader::setInputstream(std::shared_ptr<std::istream> stream) {

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

	// Load the network if the stream is available
	if (networkStream != NULL) {
		// Load the stream
		reader.setInputStream(networkStream);
		// Loop over each line of the file, which should each be PSIClusters
		do {

		} while(true/*FIXME - Conditional*/);
	}

}

