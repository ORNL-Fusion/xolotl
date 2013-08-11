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

private:

	/**
	 * The map of single-species clusters, indexed by a map that contains the
	 * name of the reactant and its size.
	 */
	std::map<std::map<std::string,int>,Reactant> singleSpeciesMap;

	/**
	 * The map of mixed or compound species clusters, indexed by a map that
	 * contains the name of the constituents of the compound reactant and their
	 * sizes.
	 */
	std::map<std::map<std::string,int>,Reactant> mixedSpeciesMap;

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

	/**
	 * This operation returns a reactant with the given name and size if it
	 * exists in the network or null if not.
	 * @param name the name of the reactant
	 * @param size the size of the reactant
	 * @return A shared pointer to the reactant
	 */
	const std::shared_ptr<Reactant> & get(std::string name, int size);

	/**
	 * This operation returns a compound reactant with the given name and size if it
	 * exists in the network or null if not.
	 * @param name the name of the compound reactant
	 * @param sizes an array containing the sizes of each piece of the reactant
	 * @return A shared pointer to the compound reactant
	 */
	const std::shared_ptr<Reactant> & getCompound(std::string name, std::vector<int> sizes);

	/**
	 * This operation adds a reactant or a compound reactant to the network.
	 * @param reactant The reactant that should be added to the network.
	 */
	void add(const std::shared_ptr<Reactant> & reactant);

	/**
	 * This operation returns the names of the reactants in the network.
	 * @return A vector with one each for each of the distinct reactant types
	 * in the network.
	 */
	const std::vector<std::string> & getNames();

	/**
	 * This operation returns the names of the compound reactants in the
	 * network.
	 * @return A vector with one each for each of the distinct compound
	 * reactant types in the network.
	 */
	const std::vector<std::string> & getCompoundNames();

	/**
	 * This operation returns a map of the properties of this reaction network.
	 * @return The map of properties that has been configured for this
	 * ReactionNetwork.
	 */
	const std::map<std::string,std::string> & getProperties();

	/**
	 * This operation sets a property with the given key to the specified value
	 * for the network. ReactionNetworks may reserve the right to ignore this
	 * operation for special key types.
	 * @param key The key for the property
	 * @param value The value to which the key should be set.
	 */
	void setProperty(std::string key, std::string value);

};

}

#endif
