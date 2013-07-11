#ifndef HEVCLUSTER_H
#define HEVCLUSTER_H

// Includes
#include "PSICluster.h"
#include <string>
#include <map>

namespace xolotlCore {

/**
 *  A cluster composed of helium and vacancies
 */
class HeVCluster : public PSICluster {

private:

	//! The number of helium atoms in this cluster.
	int numHe;

	//! The number of atomic vacancies in this cluster.
	int numV;

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	HeVCluster() : PSICluster(1) { numHe = 1; numV = 1; }

public:

	/**
	 * The constructor. All HeVClusters must be initialized with a map
	 * that describes the species of which the cluster is composed. The map
	 * should contain as its keys the names of the species and the sizes of the
	 * species as its values. The names of the species must be one of
	 * {He,V}.
	 */
	HeVCluster(int numHe, int numV);

	//! Destructor
	~HeVCluster();

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
	 * This operation returns the number of a given "species" within this
	 * cluster by passing one of {He,V} as an input argument.
	 */
	int getSpeciesSize(const std::string speciesName);

	/**
	 * This operation overrides getConnectivity() from the base class to
	 * provide the proper connectivity for a HeV cluster.
	 * @return The connectivity array
	 */
    virtual std::vector<int> getConnectivity();

    /**
	 * This operation returns the total change in this cluster due to
	 * dissociation.
	 * @return The flux due to dissociation.
	 */
	virtual double getDissociationFlux(const double temperature);

	/**
	 * This operation returns the total change in this cluster due to
	 * production.
	 * @return The flux due to this cluster being produced.
	 */
	virtual double getProductionFlux(const double temperature);

};
//end class HeVCluster

} /* end namespace xolotlCore */
#endif
