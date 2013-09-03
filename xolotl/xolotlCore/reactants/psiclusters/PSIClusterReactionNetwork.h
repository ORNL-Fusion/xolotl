#ifndef PSI_CLUSTER_REACTION_NETWORK_H
#define PSI_CLUSTER_REACTION_NETWORK_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <ReactionNetwork.h>
#include <PSICluster.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace xolotlCore {

/**
 *  This is a simple convenience class that contains a set of Reactants and the
 *  properties that describe that set.
 */
class PSIClusterReactionNetwork: public ReactionNetwork {

private:

	/**
	 * This structure compares two PSIClusters that are single species.
	 * Clusters are binned more or less uniformly between 1, 2, 3 and 4 for He,
	 * V and I respectively.
	 */
	struct PSIClusterComparator {
		bool operator()(const std::map<std::string, int>& lhs,
				const std::map<std::string, int>& rhs) const {

			// Local Declarations
			double numHe_lhs = 0, numV_lhs = 0, numI_lhs = 0;
			double numHe_rhs = 0, numV_rhs = 0, numI_rhs = 0;
			double index_lhs = 0.0, index_rhs = 0.0;
			double bigNumber = 1.0e9;

			// Get the cluster sizes, left node first
			numHe_lhs = (double) lhs.at("He");
			numV_lhs = (double) lhs.at("V");
			numI_lhs = (double) lhs.at("I");
			// Right node second
			numHe_rhs = (double) rhs.at("He");
			numV_rhs = (double) rhs.at("V");
			numI_rhs = (double) rhs.at("I");

			// Compute the indices/hashes. This simply bins the amount of each
			// time in such a way that they can be compared without a large
			// number of branches. He lands between 1 and 2, V between 2 and 3
			// and I between 3 and 4. The "big number" was chosen to be
			// sufficiently large that a single species cluster would never
			// reach that size because of physical limits.
			index_lhs = (numHe_lhs > 0) * (1.0 + (numHe_lhs / bigNumber))
					+ (numV_lhs > 0) * (2.0 + (numV_lhs / bigNumber))
					+ (numI_lhs > 0) * (3.0 + (numI_lhs / bigNumber));
			index_rhs = (numHe_rhs > 0) * (1.0 + (numHe_rhs / bigNumber))
					+ (numV_rhs > 0) * (2.0 + (numV_rhs / bigNumber))
					+ (numI_rhs > 0) * (3.0 + (numI_rhs / bigNumber));

			return index_lhs < index_rhs;
		}
	};

	/**
	 * This structure compares two PSIClusters that are mixed species. It uses
	 * a spatial hash, described in detail in the paper at:
	 * http://www.beosil.com/download/CollisionDetectionHashing_VMV03.pdf
	 *
	 * This paper can also be found by reference:
	 * Matthias Teschner et. al., "Optimized Spatial Hashing for Collision
	 * Detection of Deformable Objects," VMV 2003. Munich, Germany. November
	 * 19-21, 2003.
	 *
	 * The hash falls under the set of universal hashes for 3d vectors.
	 *
	 * The hash table size is set to 1000000 since there should never be a
	 * cluster of size one million and it is still large enough to give good
	 * hashes. Note that this number is 100x bigger than what they tested in
	 * the paper. ;)
	 */
	struct PSIMixedClusterComparator {
		bool operator()(const std::map<std::string, int>& lhs,
				const std::map<std::string, int>& rhs) const {

			// Local Declarations
			int numHe_lhs = 0, numV_lhs = 0, numI_lhs = 0;
			int numHe_rhs = 0, numV_rhs = 0, numI_rhs = 0;
			int p1 = 73856093, p2 = 19349663, p3 = 83492791;
			int hashTableSize = 10000000, hash1 = 0, hash2 = 0;

			// Get the cluster sizes
			numHe_lhs = lhs.at("He");
			numV_lhs = lhs.at("V");
			numI_lhs = lhs.at("I");
			numHe_rhs = rhs.at("He");
			numV_rhs = rhs.at("V");
			numI_rhs = rhs.at("I");

			// Compute the hashes
			hash1 =  ((numHe_lhs*p1)^(numV_lhs*p2)^(numI_lhs*p3))%hashTableSize;
			hash2 = ((numHe_rhs*p1)^(numV_rhs*p2)^(numI_rhs*p3))%hashTableSize;

//			std::cout << numHe_lhs << " " << numV_lhs << " " << numI_lhs << std::endl;
//			std::cout << numHe_rhs << " " << numV_rhs << " " << numI_rhs << std::endl;
//			std::cout << hash1 << " | " << hash2 << std::endl;

			return hash1 < hash2;
		}
	};

	/**
	 * The map of single-species clusters, indexed by a map that contains the
	 * name of the reactant and its size.
	 */
	std::map<std::map<std::string, int>, std::shared_ptr<PSICluster>,
			PSIClusterComparator> singleSpeciesMap;

	/**
	 * The map of mixed or compound species clusters, indexed by a map that
	 * contains the name of the constituents of the compound reactant and their
	 * sizes.
	 */
	std::map<std::map<std::string, int>, std::shared_ptr<PSICluster>,
			PSIMixedClusterComparator> mixedSpeciesMap;

	/**
	 * This map stores all of the clusters in the network by type
	 */
	std::map<std::string,std::shared_ptr<std::vector<std::shared_ptr<Reactant>>> > clusterTypeMap;

	/**
	 * The map of compositions to cluster ids
	 */
	std::map<std::map<std::string, int>, int> idMap;

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
	 * This operation sets the default values of the properties table and names
	 * for this network. It is used on construction and during a copy.
	 */
	void setDefaultPropsAndNames();

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
	 * This operation returns a reactant with the given name and size if it
	 * exists in the network or null if not.
	 * @param rName the name of the reactant
	 * @param size the size of the reactant
	 * @return A shared pointer to the reactant
	 */
	std::shared_ptr<Reactant> get(const std::string rName,
			const int size) const;

	/**
	 * This operation returns a compound reactant with the given name and size
	 * if it exists in the network or null if not.
	 * @param rName the name of the compound reactant
	 * @param sizes an array containing the sizes of each piece of the reactant.
	 * For PSIClusters, this array must be ordered in size by He, V and I. This
	 * array must contain an entry for He, V and I, even if only He and V or He
	 * and I are contained in the mixed-species cluster.
	 * @return A shared pointer to the compound reactant
	 */
	std::shared_ptr<Reactant> getCompound(const std::string rName,
			const std::vector<int> sizes) const;

	/**
	 * This operation returns all reactants in the network without regard for
	 * their composition or whether they are compound reactants. The list may
	 * or may not be ordered and the decision is left to implementers.
	 * @return The list of all of the reactants in the network
	 */
	std::shared_ptr<std::vector<std::shared_ptr<Reactant> > > getAll() const;

	/**
	 * This operation returns all reactants in the network with the given name.
	 * The list may or may not be ordered and the decision is left to
	 * implementers.
	 * @param name The reactant or compound reactant name
	 * @return The list of all of the reactants in the network or null if the
	 * name is invalid.
	 */
	std::shared_ptr<std::vector<std::shared_ptr<Reactant> > > getAll(std::string name) const;

	/**
	 * This operation adds a reactant or a compound reactant to the network.
	 * Adding a reactant to the network does not set the network as the
	 * reaction network for the reactant. This step must be performed
	 * separately to allow for the scenario where the network is generated
	 * entirely before running.
	 *
	 * The reactant will not be added to the network if the PSICluster does
	 * not recognize it as a type of reactant that it cares about (including
	 * adding null).
	 * @param reactant The reactant that should be added to the network.
	 */
	void add(std::shared_ptr<Reactant> reactant);

	/**
	 * This operation returns the names of the reactants in the network. For a
	 * PSIClusterReactionNetwork, these are He, V, I, HeV, HeI.
	 * @return A vector with one entry for each of the distinct reactant types
	 * in the network.
	 */
	const std::vector<std::string> & getNames() const;

	/**
	 * This operation returns the names of the compound reactants in the
	 * network.
	 * @return A vector with one each for each of the distinct compound
	 * reactant types in the network.
	 */
	const std::vector<std::string> & getCompoundNames() const;

	/**
	 * This operation returns a map of the properties of this reaction network.
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
	 * @param key The key for the property
	 * @param value The value to which the key should be set.
	 */
	void setProperty(std::string key, std::string value);

	/**
	 * This operation returns the size or number of reactants in the network.
	 * @return The number of reactants in the network
	 */
	int size();

	/**
	 * This operation returns the id of a reactant if it exists in the network.
	 * @param reactant The reactant
	 * @return The id of the reactant. This id is guaranteed to be between 1 and
	 * n, including both, for n reactants in the network.
	 */
	int getReactantId(const Reactant & reactant);

	/**
	 * This is a utility operation that creates a composition vector with an
	 * entry for helium, vacancies and interstitials. It is handy because it
	 * removes the need to construct the vectors properly and locally. The
	 * vector can be used to retrieve clusters.
	 *
	 * This function will never return an composition with less than three
	 * elements and it will always return element sizes greater than zero. If
	 * any of the elements are negative, it will default to 1,0,0 (single He).
	 * @param numHe The number of helium atoms in the cluster
	 * @param numV The number of atomic vacancies in the cluster
	 * @param numI The number of interstitial defects in the cluster
	 * @return A vector of size three with an entry for each of the parts
	 * equal to the numbers that were passed for that part.
	 */
	std::vector<int> getCompositionVector(int numHe, int numV, int numI) {
		// This flag is used so that negative numbers can be checked without
		// branching.
		int hasNegativeElement = ((numHe < 0) + (numV < 0) + (numI < 0) > 0);
		std::vector<int> composition(3);
		composition[0] = std::max(1*hasNegativeElement,numHe);
		composition[1] = std::max(0,numV);
		composition[2] = std::max(0,numI);
		return composition;
	}

};

}

#endif
