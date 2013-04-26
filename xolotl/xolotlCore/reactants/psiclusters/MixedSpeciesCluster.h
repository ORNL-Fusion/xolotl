#ifndef MIXEDSPECIESCLUSTER_H
#define MIXEDSPECIESCLUSTER_H

// Includes
#include "PSICluster.h"
#include <string>
#include <map>

namespace xolotlCore {

/**
 *  This class represents a cluster composed of multiple atomic species,
 *  vacancies and interstitials.
 */
class MixedSpeciesCluster : public PSICluster {

private:

	//! The number of hydrogen atoms in this cluster.
	int numH;

	//! The number of interstitial defects in this cluster.
	int numI;

	//! The number of deuterium atoms in this cluster.
	int numD;

	//! The number of helium atoms in this cluster.
	int numHe;

	//! The number of tritium atoms in this cluster.
	int numT;

	//! The number of atomic vacancies in this cluster.
	int numV;

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	MixedSpeciesCluster():PSICluster(1) {}

public:

	/**
	 * The constructor. All MixedSpeciesClusters must be initialized with a map
	 * that describes the species of which the cluster is composed. The map
	 * should contain as its keys the names of the species and the sizes of the
	 * species as its values. The names of the species must be one of
	 * {H,He,I,V,D,T}.
	 */
	MixedSpeciesCluster(const std::map<std::string,int> speciesMap);

	//! Destructor
	~MixedSpeciesCluster();

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
	 * cluster by passing one of {H,He,I,V,D,T} as an input argument.
	 */
	int getSpeciesSize(const std::string speciesName);

};
//end class MixedSpeciesCluster

} /* end namespace xolotlCore */
#endif
