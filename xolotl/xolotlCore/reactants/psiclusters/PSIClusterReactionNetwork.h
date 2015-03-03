#ifndef PSI_CLUSTER_REACTION_NETWORK_H
#define PSI_CLUSTER_REACTION_NETWORK_H

// Includes
//#include <xolotlPerf.h>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <ReactionNetwork.h>
#include <PSICluster.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <Constants.h>

// Override the hash operation for the composition maps used by the
// PSIClusterReactionNetwork to store reactants.
namespace std {
template<>
class hash<std::map<std::string, int>> {
public:
	long operator()(const std::map<std::string, int>& composition) const {
		int bigNumber = 1e9;
		return (composition.at(xolotlCore::heType) * 10 + composition.at(xolotlCore::vType) * 200
				+ composition.at(xolotlCore::iType) * 3000) * bigNumber;
	}
};
}

namespace xolotlCore {

/**
 *  This class manages the set of reactants and compound reactants (
 *  combinations of normal reactants) for PSI clusters. It also manages a
 *  set of properties that describes the total collection.
 *
 *  This class is a very heavyweight class that should not be abused.
 *
 *  Reactants that are added to this network must be added as with shared_ptrs.
 *  Furthermore, reactants that are added to this network have their ids set to
 *  a network specific id. Reactants should not be shared between separate
 *  instances of a PSIClusterReactionNetwork.
 */
class PSIClusterReactionNetwork: public ReactionNetwork {

private:

	/**
	 * The map of single-species clusters, indexed by a map that contains the
	 * name of the reactant and its size.
	 */
	std::unordered_map< std::map< std::string, int >,
		std::shared_ptr<PSICluster> > singleSpeciesMap;

	/**
	 * The map of mixed or compound species clusters, indexed by a map that
	 * contains the name of the constituents of the compound reactant and their
	 * sizes.
	 */
	std::unordered_map< std::map< std::string, int >,
		std::shared_ptr<PSICluster> > mixedSpeciesMap;

	/**
	 * This map stores all of the clusters in the network by type.
	 */
	std::map<std::string, std::shared_ptr<
		std::vector< std::shared_ptr<Reactant> > > > clusterTypeMap;

	/**
	 * The names of the reactants supported by this network.
	 */
	std::vector<std::string> names;

	/**
	 * The names of the compound reactants supported by this network.
	 */
	std::vector<std::string> compoundNames;

	/**
	 * The size of the network. It is also used to set the id of new Reactants
	 * that are added to the network.
	 */
	int networkSize;

	/**
	 * The current temperature at which the network's clusters exist.
	 */
	double temperature;

	/**
	 * The list of all of the reactants in the network. This list is filled and
	 * maintained by the getAll() operation.
	 */
	std::shared_ptr< std::vector<Reactant *> > allReactants;

	/**
	 * This operation sets the default values of the properties table and names
	 * for this network. It is used on construction and during a copy.
	 */
	void setDefaultPropsAndNames();

	/**
	 * The Constructor
	 */
	PSIClusterReactionNetwork();

public:

	/**
	 * The Constructor
	 *
	 * @param registry The performance handler registry
	 */
	PSIClusterReactionNetwork(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * The copy constructor.
	 *
	 * @param other
	 */
	PSIClusterReactionNetwork(const PSIClusterReactionNetwork &other);

	/**
	 * This operation sets the temperature at which the reactants currently
	 * exists. It calls setTemperature() on each reactant.
	 *
	 * This is the simplest way to set the temperature for all reactants is to
	 * call the ReactionNetwork::setTemperature() operation.
	 *
	 * @param temp The new temperature
	 */
	virtual void setTemperature(double temp);

	/**
	 * This operation returns the temperature at which the cluster currently exists.
	 *
	 * @return The temperature.
	 */
	virtual double getTemperature() const;

	/**
	 * This operation returns a reactant with the given type and size if it
	 * exists in the network or null if not.
	 *
	 * @param type the type of the reactant
	 * @param size the size of the reactant
	 * @return A pointer to the reactant
	 */
	Reactant * get(const std::string type,
			const int size) const;

	/**
	 * This operation returns a compound reactant with the given type and size
	 * if it exists in the network or null if not.
	 *
	 * @param type the type of the compound reactant
	 * @param sizes an array containing the sizes of each piece of the reactant.
	 * For PSIClusters, this array must be ordered in size by He, V and I. This
	 * array must contain an entry for He, V and I, even if only He and V or He
	 * and I are contained in the mixed-species cluster.
	 * @return A pointer to the compound reactant
	 */
	Reactant * getCompound(const std::string type,
			const std::vector<int> sizes) const;

	/**
	 * This operation returns all reactants in the network without regard for
	 * their composition or whether they are compound reactants. The list may
	 * or may not be ordered and the decision is left to implementers.
	 *
	 * @return The list of all of the reactants in the network
	 */
	const std::shared_ptr<std::vector<Reactant *>> & getAll() const;

	/**
	 * This operation returns all reactants in the network with the given name.
	 * The list may or may not be ordered and the decision is left to
	 * implementers.
	 *
	 * @param name The reactant or compound reactant name
	 * @return The list of all of the reactants in the network or null if the
	 * name is invalid.
	 */
	std::vector<Reactant *> getAll(std::string name) const;

	/**
	 * This operation adds a reactant or a compound reactant to the network.
	 * Adding a reactant to the network does not set the network as the
	 * reaction network for the reactant. This step must be performed
	 * separately to allow for the scenario where the network is generated
	 * entirely before running.
	 *
	 * This operation sets the id of the reactant to one that is specific
	 * to this network. Do not share Reactants across networks! This id is
	 * guaranteed to be between 1 and n, including both, for n reactants in
	 * the network.
	 *
	 * The reactant will not be added to the network if the PSICluster does
	 * not recognize it as a type of reactant that it cares about (including
	 * adding null). This operation throws an exception of type std::string
	 * if the reactant is  already in the network.
	 *
	 * @param reactant The reactant that should be added to the network.
	 */
	void add(std::shared_ptr<Reactant> reactant);

	/**
	 * This operation returns the names of the reactants in the network. For a
	 * PSIClusterReactionNetwork, these are He, V, I, HeV, HeI.
	 *
	 * @return A vector with one entry for each of the distinct reactant types
	 * in the network
	 */
	const std::vector<std::string> & getNames() const;

	/**
	 * This operation returns the names of the compound reactants in the
	 * network.
	 *
	 * @return A vector with one each for each of the distinct compound
	 * reactant types in the network
	 */
	const std::vector<std::string> & getCompoundNames() const;

	/**
	 * This operation returns a map of the properties of this reaction network.
	 *
	 * @return The map of properties that has been configured for this
	 * ReactionNetwork.
	 *
	 * The PSIClusterReactionNetwork always has the following properties:
	 * > maxHeClusterSize - The number of He atoms in the largest single-species
	 *  He cluster.
	 * > maxVClusterSize - The number of atomic vacancies in the largest
	 * single-species V cluster.
	 * > maxIClusterSize - The number of interstitials in the largest
	 * single-species I cluster.
	 * > maxHeVClusterSize - The number of species of all types in the largest
	 * mixed species in the network. It is equal to the sum of the max single
	 * species helium and vacancy cluster sizes by default.
	 * > maxHeIClusterSize - The number of species of all types in the largest
	 * mixed species in the network. It is equal to the sum of the max single
	 * species helium and vacancy cluster sizes by default.
	 * > numHeClusters - The number of single-species He clusters of all sizes in
	 * the network.
	 * > numVClusters - The number of single-species V clusters of all sizes in the
	 * network.
	 * > numIClusters - The number of single-species I clusters of all sizes in the
	 * network.
	 * > numHeVClusters - The number of HeV clusters of all sizes in the
	 * network.
	 * > numHeIClusters - The number of HeI clusters of all sizes in the
	 * network.
	 *
	 * These properties are always updated when a cluster is added.
	 */
	const std::map<std::string, std::string> & getProperties();

	/**
	 * This operation sets a property with the given key to the specified value
	 * for the network. ReactionNetworks may reserve the right to ignore this
	 * operation for special key types, most especially those that they manage
	 * on their own.
	 *
	 * @param key The key for the property
	 * @param value The value to which the key should be set
	 */
	void setProperty(std::string key, std::string value);

	/**
	 * This operation returns the size or number of reactants in the network.
	 *
	 * @return The number of reactants in the network
	 */
	int size();

	/**
	 * This is a utility operation that creates a composition vector with an
	 * entry for helium, vacancies and interstitials. It is handy because it
	 * removes the need to construct the vectors properly and locally. The
	 * vector can be used to retrieve clusters.
	 *
	 * This function will never return a composition with less than three
	 * elements and it will always return element sizes greater than zero. If
	 * any of the elements are negative, it will default to 1,0,0 (single He).
	 *
	 * @param numHe The number of helium atoms in the cluster
	 * @param numV The number of atomic vacancies in the cluster
	 * @param numI The number of interstitial defects in the cluster
	 *
	 * @return A vector of size three with an entry for each of the parts
	 * equal to the numbers that were passed for that part.
	 */
	std::vector<int> getCompositionVector(int numHe, int numV, int numI) {
		// This flag is used so that negative numbers can be checked without
		// branching.
		int hasNegativeElement = ((numHe < 0) + (numV < 0) + (numI < 0) > 0);

		std::vector<int> composition(3);

		composition[0] = std::max(0,numHe*!hasNegativeElement);
		composition[1] = std::max(0,numV*!hasNegativeElement);
		composition[2] = std::max(0,numI*!hasNegativeElement);

		return composition;
	}
};

}

#endif
