#ifndef SUPERCLUSTER_H
#define SUPERCLUSTER_H

// Includes
#include "PSICluster.h"
#include <string>

namespace xolotlCore {

/**
 *  A cluster gathering the average properties of many HeV clusters.
 */
class SuperCluster: public PSICluster {

private:

	//! The mean number of helium atoms in this cluster.
	double numHe;

	//! The mean number of atomic vacancies in this cluster.
	double numV;

	//! The number of HeV clusters gathered in this one.
	int sectionWidth;

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	SuperCluster() :
		PSICluster(1)
	{ numHe = 1.0; numV = 1.0; }

public:

	//! The vector of HeV clusters it will replace
	std::vector<PSICluster *> heVVector;

	/**
	 * The constructor. All SuperClusters must be initialized with its
	 * composition.
	 *
	 * @param numHe The mean number of helium atoms in this cluster
	 * @param numV The mean number of vacancies in this cluster
	 * @param width The width of this super cluster
	 * @param radius The mean radius
	 * @param energy The mean formation energy
	 * @param registry The performance handler registry
	 */
	SuperCluster(double numHe, double numV, int width, double radius, double energy,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor.
	 *
	 * @param other the reactant to be copied
	 */
	SuperCluster(const SuperCluster &other);

	//! Destructor
	~SuperCluster() {}

	/**
	 * This operation returns a Reactant that is created using the copy
	 * constructor of SuperCluster.
	 *
	 * @return A copy of this reactant
	 */
	virtual std::shared_ptr<Reactant> clone();

	/**
	 * This operation returns true to signify that this cluster is a mixture of
	 * He and V.
	 *
	 * @return True if mixed
	 */
	virtual bool isMixed() const {return true;}

	/**
	 * Set the HeV vector
	 */
	void setHeVVector(std::vector<PSICluster *> vec)
		{heVVector = vec;}

	/**
	 * Computes a row of the reaction connectivity matrix corresponding to
	 * this reactant.
	 *
	 * If two reactants alone can form a reaction, the element at the position
	 * of the second reactant is 1, otherwise 0.
	 */
	void createReactionConnectivity();

	/**
	 * Computes a row of the dissociation connectivity matrix corresponding to
	 * this reactant.
	 *
	 * If two reactants together can be produced by a single reaction,
	 * the element at the position of the second reactant is 1, otherwise 0.
	 */
	void createDissociationConnectivity();

	/**
	 * Calculate all the rate constants for the reactions and dissociations in which this
	 * cluster is taking part. Store these values in the kConstant field of ClusterPair
	 * or CombiningCluster. Need to be called only when the temperature changes.
	 */
	void computeRateConstants();

};
//end class SuperCluster

} /* end namespace xolotlCore */
#endif
