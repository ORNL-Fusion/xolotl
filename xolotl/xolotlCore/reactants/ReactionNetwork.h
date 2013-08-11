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
 *  This class manages the set of reactant and compound reactants (
 *  combinations of normal reactants). It also manages a set of properties
 *  that describe both.
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
	std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>reactants;

	/** The constructor. It initializes the properties map and reactants vector.
	 */
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
	virtual std::map<std::string, int> toClusterMap(int index) const;

	/**
	 * Converts an cluster map (with `speciesLabel` => `quantity`)
	 * to the index corresponding to its position in the reactants vector.
	 * This operation is const, but the cluster map can not be marked const
	 * because of the mechanics of effectively using maps in C++. Sorry
	 * idealists! ;)
	 */
	virtual int toClusterIndex(std::map<std::string, int> clusterMap) const;

	/**
	 * This operation returns a reactant with the given name and size if it
	 * exists in the network or null if not.
	 * @param name the name of the reactant
	 * @param size the size of the reactant
	 * @return A shared pointer to the reactant
	 */
	virtual const std::shared_ptr<Reactant> & get(std::string name, int size) = 0;

	/**
	 * This operation returns a compound reactant with the given name and size if it
	 * exists in the network or null if not.
	 * @param name the name of the compound reactant
	 * @param sizes an array containing the sizes of each piece of the reactant
	 * @return A shared pointer to the compound reactant
	 */
	virtual const std::shared_ptr<Reactant> & getCompound(std::string name, std::vector<int> sizes) = 0;

	/**
	 * This operation adds a reactant or a compound reactant to the network.
	 * @param reactant The reactant that should be added to the network.
	 */
	virtual void add(const std::shared_ptr<Reactant> & reactant) = 0;

	/**
	 * This operation returns the names of the reactants in the network.
	 * @return A vector with one each for each of the distinct reactant types
	 * in the network.
	 */
	virtual const std::vector<std::string> & getNames() = 0;

	/**
	 * This operation returns the names of the compound reactants in the
	 * network.
	 * @return A vector with one each for each of the distinct compound
	 * reactant types in the network.
	 */
	virtual const std::vector<std::string> & getCompoundNames() = 0;

	/**
	 * This operation returns a map of the properties of this reaction network.
	 * @return The map of properties that has been configured for this
	 * ReactionNetwork.
	 */
	virtual const std::map<std::string,std::string> & getProperties() = 0;

	/**
	 * This operation sets a property with the given key to the specified value
	 * for the network. ReactionNetworks may reserve the right to ignore this
	 * operation for special key types.
	 * @param key The key for the property
	 * @param value The value to which the key should be set.
	 */
	virtual void setProperty(std::string key, std::string value) = 0;

	/**
	 * The destructor
	 */
	virtual ~ReactionNetwork() {}

};

}

#endif
