/*
 * PSIClusterNetworkLoader.cpp
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#include "psiclusternetworkloader.h"

namespace xolotlCore {

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
}

/**
 * This operation will load the reaction network from the inputstream in
 * the format specified previously. The network will be empty if it can not
 * be loaded.
 * @return The reaction network //Should it be vector<vector>?
 */
std::vector<PSICluster> PSIClusterNetworkLoader::load() {

	std::vector<PSICluster> vector;

	return vector;
}

/**
 * This operation returns the properties map that was created for the
 * network. If it is called before the network is loaded or the network can
 * not be loaded, it will be empty.
 * @return The property map
 */
std::map<std::string, std::string> PSIClusterNetworkLoader::getProperties() {

	std::map<std::string, std::string> map;

	return map;
}

} /* namespace xolotlCore */
