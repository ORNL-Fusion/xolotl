#ifndef HEVCLUSTER_H
#define HEVCLUSTER_H

// Includes
#include "PSICluster.h"
#include "../../xolotlPerf/HandlerRegistryFactory.h"
#include <string>
#include <map>

namespace xolotlCore {

/**
 *  A cluster composed of helium and vacancies
 */
class HeVCluster: public PSICluster {

private:

	//! The helium cluster of size 1
	PSICluster * heCluster;

	//! The vacancy cluster of size 1
	PSICluster * vCluster;

	//! The interstitial cluster of size 1
	PSICluster * iCluster;

	//! The HeV cluster with one less helium
	PSICluster * heVClusterLessHe;

	//! The HeV cluster with one less vacancy
	PSICluster * heVClusterLessV;

	//! The number of helium atoms in this cluster.
	int numHe;

	//! The number of atomic vacancies in this cluster.
	int numV;

	//! The sum of the dissociation constants between this cluster and the
	//! clusters of size 1.
	double f4 = 0.0;

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	HeVCluster() :
		PSICluster(1)
	{ numHe = 1; numV = 1; }

protected:

	/**
	 * This operation computes the partial derivatives due to dissociation
	 * reactions. The partial derivatives due to dissociation for compound
	 * clusters are significantly different than those single-species clusters.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param temperature The temperature at which the reactions are occurring.
	 */
	virtual void getDissociationPartialDerivatives(std::vector<double> & partials,
			double temperature) const;

public:

	/**
	 * The constructor. All HeVClusters must be initialized with a map
	 * that describes the species of which the cluster is composed. The map
	 * should contain as its keys the names of the species and the sizes of the
	 * species as its values. The names of the species must be one of
	 * {He,V}.
	 */
	HeVCluster(int numHe, int numV,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor
	 * @param other the reactant to be copied
	 */
	HeVCluster(const HeVCluster &other);

	//! Destructor
	~HeVCluster();

	/**
	 * This operation returns a Reactant that is created using the copy
	 * constructor of HeVCluster.
	 * @return A copy of this reactant.
	 */
	virtual std::shared_ptr<Reactant> clone();

	/**
	 * This operation returns the total generation rate due to emission for
	 * this cluster.
	 */
	double getGenByEm();

	/**
	 * This operation returns the total annihilation rate due to emission for
	 * this cluster.
	 */
	double getAnnByEm();

	/**
	 * This operation returns the total change in this cluster due to
	 * dissociation.
	 * @return The flux due to dissociation.
	 */
	virtual double getDissociationFlux(double temperature) const;

	/**
	 * This operation returns true to signify that this cluster is a mixture of
	 * He and V.
	 */
	virtual bool isMixed() const { return true; };

	/**
	 * This operation overrides Reactant's setTemperature operation to
	 * correctly recompute the diffusion coefficient and other
	 * temperature-dependent quantities when the temperature is set.
	 * @param temp
	 */
	virtual void setTemperature(double temp);

protected:

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

};
//end class HeVCluster

} /* end namespace xolotlCore */
#endif
