#ifndef REACTION_NETWORK_H
#define REACTION_NETWORK_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include <map>

namespace xolotlCore {

class Reactant;

/**
 *  This is a simple convenience class that contains a set of Reactants and the
 *  properties that describe that set.
 */
class ReactionNetwork {
public:

	/**
	 * The properties of this network. The exact configuration of the map is
	 * specified by the class that loaded the network.
	 */
	std::shared_ptr<std::map<std::string, std::string>> properties;

	/**
	 * The set of reactants. The exact order of the reactants in the set is
	 * specified by the class that loaded them.
	 */
	std::shared_ptr<std::vector<std::shared_ptr<Reactant>>> reactants;

	//! The constructor. It initializes the properties map and reactants vector.
	ReactionNetwork() :
	properties(new std::map<std::string, std::string>()),
	reactants(new std::vector<std::shared_ptr<Reactant>>())
	{}
	
	/**
	 * The copy constructor
	 * @param other The ReactionNetwork to copy
	 */
	ReactionNetwork(const ReactionNetwork &other);
	
	/**
	 * Converts an cluster index (found in the `reactants` vector)
	 * to a map describing the cluster's 
	 *
	 * @returns a map with `speciesLabel` => `quantity`
	 */
	virtual std::map<std::string, int> toClusterMap(int index);
	
	/**
	 * Converts an cluster map (with `speciesLabel` => `quantity`)
	 * to the index corresponding to its position in the reactants vector
	 */
	virtual int toClusterIndex(std::map<std::string, int> clusterMap);
	
	/**
	 * The destructor
	 */
	virtual ~ReactionNetwork() {}

	/**
	 * Return whether reactantI reacts with (or is connected in the network)
	 * reactantJ.
	 * @param reactantI
	 * @param reactantJ
	 * @return
	 */
	bool isConnected(int reactantI, int reactantJ);

};

}

#endif
