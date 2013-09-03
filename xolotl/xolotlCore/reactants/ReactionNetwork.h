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

protected:

	/**
	 * The properties of this network. The exact configuration of the map is
	 * specified by the class that loaded the network.
	 */
	std::shared_ptr<std::map<std::string, std::string>> properties;

public:

	/** The constructor. It initializes the properties map and reactants vector.
	 */
	ReactionNetwork() :
	properties(new std::map<std::string, std::string>())
	{}

	/**
	 * The copy constructor
	 * @param other The ReactionNetwork to copy
	 */
	ReactionNetwork(const ReactionNetwork &other);

	/**
	 * The destructor
	 */
	virtual ~ReactionNetwork() {}

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
	virtual std::shared_ptr<Reactant> get(const std::string rName, const int size) const = 0;

	/**
	 * This operation returns a compound reactant with the given name and size if it
	 * exists in the network or null if not.
	 * @param name the name of the compound reactant
	 * @param sizes an array containing the sizes of each piece of the reactant
	 * @return A shared pointer to the compound reactant
	 */
	virtual std::shared_ptr<Reactant> getCompound(const std::string rName, const std::vector<int> sizes) const = 0;

	/**
	 * This operation returns all reactants in the network without regard for
	 * their composition or whether they are compound reactants. The list may
	 * or may not be ordered and the decision is left to implementers.
	 * @return The list of all of the reactants in the network
	 */
	virtual std::shared_ptr<std::vector<std::shared_ptr<Reactant> > > getAll() const = 0;

	/**
	 * This operation returns all reactants in the network with the given name.
	 * The list may or may not be ordered and the decision is left to
	 * implementers.
	 * @param name The reactant or compound reactant name
	 * @return The list of all of the reactants in the network or null if the
	 * name is invalid.
	 */
	virtual std::shared_ptr<std::vector<std::shared_ptr<Reactant> > > getAll(std::string name) const = 0;

	/**
	 * This operation adds a reactant or a compound reactant to the network.
	 * Adding a reactant to the network does not set the network as the
	 * reaction network for the reactant. This step must be performed
	 * separately to allow for the scenario where the network is generated
	 * entirely before running.
	 * @param reactant The reactant that should be added to the network.
	 */
	virtual void add(std::shared_ptr<Reactant> reactant) = 0;

	/**
	 * This operation returns the names of the reactants in the network.
	 * @return A vector with one entry for each of the distinct reactant types
	 * in the network.
	 */
	virtual const std::vector<std::string> & getNames() const = 0;

	/**
	 * This operation returns the names of the compound reactants in the
	 * network.
	 * @return A vector with one each for each of the distinct compound
	 * reactant types in the network.
	 */
	virtual const std::vector<std::string> & getCompoundNames() const = 0;

	/**
	 * This operation returns a map of the properties of this reaction network.
	 * @return The map of properties that has been configured for this
	 * ReactionNetwork.
	 */
	virtual const std::map<std::string,std::string> & getProperties() = 0;

	/**
	 * This operation sets a property with the given key to the specified value
	 * for the network. ReactionNetworks may reserve the right to ignore this
	 * operation for special key types, most especially those that they manage
	 * on their own.
	 * @param key The key for the property
	 * @param value The value to which the key should be set.
	 */
	virtual void setProperty(std::string key, std::string value) = 0;

	/**
	 * This operation returns the size or number of reactants in the network.
	 * @return The number of reactants in the network
	 */
	virtual int size() = 0;

	/**
	 * This operation returns the id of a reactant if it exists in the network.
	 * @param reactant The reactant
	 * @return The id of the reactant. This id is guaranteed to be between 1 and
	 * n, including both, for n reactants in the network.
	 */
	virtual int getReactantId(const Reactant & reactant) = 0;

};

}

#endif
