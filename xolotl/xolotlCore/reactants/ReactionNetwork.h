#ifndef REACTION_NETWORK_H
#define REACTION_NETWORK_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include <map>

namespace xolotlPerf {
	class IHandlerRegistry;
    class IEventCounter;
};

namespace xolotlCore {

class Reactant;

/**
 *  This class manages the set of reactants and compound reactants (
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

    /**
     * The performance handler registry that will be used with
     * this class.
     */
    std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry;

    /**
     * Counter for the number of times the network concentration is updated.
     */
    std::shared_ptr<xolotlPerf::IEventCounter> concUpdateCounter;

	/**
	 * The default constructor. It initializes the properties map and reactants vector.
	 */
	ReactionNetwork();

public:

	/**
	 * The constructor that takes the performance handler registry.
	 * It initializes the properties map and reactants vector.
	 */
	ReactionNetwork(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);


	/**
	 * The copy constructor
	 * @param other The ReactionNetwork to copy
	 */
	ReactionNetwork(const ReactionNetwork &other);

	/**
	 * The destructor
	 */
	virtual ~ReactionNetwork();

	/**
	 * This operation sets the temperature at which the reactants currently
	 * exists. It calls setTemperature() on each reactant.
	 *
	 * This is the simplest way to set the temperature for all reactants is to
	 * call the ReactionNetwork::setTemperature() operation.
	 *
	 * @param temp The new temperature
	 */
	virtual void setTemperature(double temp) = 0;

	/**
	 * This operation returns the temperature at which the cluster currently exists.
	 * @return The temperature.
	 */
	virtual double getTemperature() const = 0;

	/**
	 * This operation returns a reactant with the given name and size if it
	 * exists in the network or null if not.
	 * @param type the type of the reactant
	 * @param size the size of the reactant
	 * @return A shared pointer to the reactant
	 */
	virtual Reactant * get(const std::string type, const int size) const = 0;

	/**
	 * This operation returns a compound reactant with the given name and size if it
	 * exists in the network or null if not.
	 * @param type the type of the compound reactant
	 * @param sizes an array containing the sizes of each piece of the reactant
	 * @return A shared pointer to the compound reactant
	 */
	virtual Reactant * getCompound(const std::string type, const std::vector<int> sizes) const = 0;

	/**
	 * This operation returns all reactants in the network without regard for
	 * their composition or whether they are compound reactants. The list may
	 * or may not be ordered and the decision is left to implementers.
	 * @return The list of all of the reactants in the network
	 */
	virtual const std::shared_ptr<std::vector<Reactant *>> & getAll() const = 0;

	/**
	 * This operation returns all reactants in the network with the given name.
	 * The list may or may not be ordered and the decision is left to
	 * implementers.
	 * @param name The reactant or compound reactant name
	 * @return The list of all of the reactants in the network or null if the
	 * name is invalid.
	 */
	virtual std::vector<Reactant *> getAll(std::string name) const = 0;

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
	 * This operation fills an array of doubles with the concentrations of all
	 * of the Reactants in the network.
	 * @param concentrations The array that will be filled with the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. If the array is to small to hold the concentrations, SIGSEGV will
	 * be thrown.
	 */
	void fillConcentrationsArray(double * concentrations);

	/**
	 * This operation updates the concentrations for all Reactants in the
	 * network from an array.
	 * @param concentrations The array of doubles that will be for the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. Properly aligning the array in memory so that this operation
	 * does not overrun is up to the caller.
	 */
	void updateConcentrationsFromArray(double * concentrations);



    /**
     * Request that all Reactants in the network release their 
     * pointers to the network, to break cycles and allow the
     * network and the clusters/reactant objects it contains to 
     * be destroyed gracefully.
     *
     * Should only be done when the network is no longer needed.
     * Ideally, we would do this from the network's destructor.
     * However, unless the reactants in the network release their
     * shared_ptrs to the network, the reference count on the
     * network's shared_ptr will never reach zero and the 
     * object it owns will never be destroyed.  (Hence, the
     * reactant objects will also never be destroyed.)
     */
    void askReactantsToReleaseNetwork( void );
};

}

#endif
