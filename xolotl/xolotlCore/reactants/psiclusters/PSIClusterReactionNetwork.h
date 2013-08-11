#ifndef PSI_CLUSTER_REACTION_NETWORK_H
#define PSI_CLUSTER_REACTION_NETWORK_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include <map>
#include "../ReactionNetwork.h"

namespace xolotlCore {

/**
 *  This is a simple convenience class that contains a set of Reactants and the
 *  properties that describe that set.
 */
class PSIClusterReactionNetwork : public ReactionNetwork {

public:

	/**
	 * The Constructor
	 */
	PSIClusterReactionNetwork();

	/**
	 * The copy constructor
	 * @param other
	 */
	PSIClusterReactionNetwork(const PSIClusterReactionNetwork &other);
	
	/**
	 * Converts an cluster index (found in the `reactants` vector)
	 * to a map describing the cluster's 
	 *
	 * @returns a map with `speciesLabel` => `quantity`
	 */
	std::map<std::string, int> toClusterMap(int index) const;
	
	/**
	 * Converts an cluster map (with `speciesLabel` => `quantity`)
	 * to the index corresponding to its position in the reactants vector
	 * @param the map of species labels to species size that describes the
	 * desired cluster, i.e. "He",1 or "He",1;"V",2.
	 * @returns the index of the cluster or -1 if it could not be found
	 */
	int toClusterIndex(std::map<std::string, int> clusterMap) const;
};

}

#endif
