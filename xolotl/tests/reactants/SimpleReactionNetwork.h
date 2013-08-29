#ifndef SIMPLEREACTIONNETWORK_H_
#define SIMPLEREACTIONNETWORK_H_

#include <PSIClusterReactionNetwork.h>

namespace testUtils {

/**
 * This class creates a simple reaction network used for testing. It contains
 * 10 each of He, V and I clusters and twenty HeV clusters. The HeV clusters are
 * stored in ascending order of He and then V, (1,1;2,1;2,2;3,2;3,3;etc.).
 *
 * It does not register itself as the ReactionNetwork for its clusters because
 * of limitations with shared_ptrs and "this." So, the
 * TestUtils::getSimpleReactionNetwork() operation should always be called to
 * insure that it is properly initialized.
 */
class SimpleReactionNetwork : public xolotlCore::PSIClusterReactionNetwork {

public:
	//! Constructor
	SimpleReactionNetwork();

	//! Destructor
	virtual ~SimpleReactionNetwork();

	/**
		 * This operation returns a reactant with the given name and size if it
		 * exists in the network or null if not.
		 * @param name the name of the reactant
		 * @param size the size of the reactant
		 * @return A shared pointer to the reactant
		 */
		std::shared_ptr<xolotlCore::Reactant> get(const std::string rName, const int size);

		/**
		 * This operation returns a compound reactant with the given name and size if it
		 * exists in the network or null if not.
		 * @param name the name of the compound reactant
		 * @param sizes an array containing the sizes of each piece of the reactant
		 * @return A shared pointer to the compound reactant
		 */
		std::shared_ptr<xolotlCore::Reactant> getCompound(const std::string rName, const std::vector<int> sizes);

		/**
		 * This operation adds a reactant or a compound reactant to the network.
		 * @param reactant The reactant that should be added to the network.
		 */
		void add(std::shared_ptr<xolotlCore::Reactant> reactant);

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
		void setProperty(const std::string key, const std::string value);

		/**
		 * This operation returns all reactants in the network without regard for
		 * their composition or whether they are compound reactants. The list may
		 * or may not be ordered and the decision is left to implementers.
		 * @return The list of all of the reactants in the network
		 */
		std::shared_ptr<std::vector<std::shared_ptr<xolotlCore::Reactant> > > getAll();

};

/**
 * This operation creates a SimpleReactionNetwork and makes sure that it is
 * properly registered with the clusters it contains. This operation should
 * always be called instead of constructing a SimpleReactionNetwork manually.
 * @return The reaction network.
 */
std::shared_ptr<xolotlCore::ReactionNetwork>  getSimpleReactionNetwork();

} /* end namespace testUtils */
#endif
