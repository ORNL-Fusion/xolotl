#ifndef PSICLUSTER_H
#define PSICLUSTER_H

// Includes
#include "Reactant.h"

namespace xolotlCore {

/**
 * The PSICluster class is a Reactant that is specialized to work for
 * simulations of plasma-surface interactions. It provides special routines
 * for calculating the total flux due to production and dissociation and
 * obtaining the cluster size.
 *
 * PSIClusters must always be initialized with a size. If the constructor is
 * passed a size of zero or less, the actual size will be set to 1.
 */
class PSICluster: public Reactant {

protected:

	/**
	 * The total size of this cluster including the contributions from all
	 * species.
	 */
	int size;

	/**
	 * The diffusion factor, D_0, that is used to calculate the diffusion
	 * coefficient for this cluster. The default value is 0 (does not diffuse).
	 */
	double diffusionFactor;

	/**
	 * The binding energies for this cluster with clusters of other species
	 * types. There is one binding energy for each of the other species ordered
	 * by He, V, I and mixed species at indices 0, 1, 2 and 3 respectively.
	 */
	double bindingEnergies[4];

public:

	//! Constructor
	PSICluster(int clusterSize);

	//! Destructor
	virtual ~PSICluster();

	/**
	 * This operation returns the total change in this cluster due to
	 * dissociation.
	 * @return The flux due to dissociation.
	 */
	virtual double getDissociationFlux();

	/**
	 * This operation returns the total change in this cluster due to
	 * production.
	 * @return The flux due to this cluster being produced.
	 */
	virtual double getProductionFlux();

	/**
	 * This operation returns the total size of the cluster.
	 * @return The total size of this cluster including the contributions
	 * from all species types.
	 */
	virtual int getSize();

	/**
	 * This operation returns the total generation rate of this cluster due
	 * to captures.
	 * @return The total amount of this cluster generated due to captures.
	 */
	virtual double getGenByCapt();

	/**
	 * This operation returns the total generation rate of this cluster due
	 * to annihilation.
	 * @return The total amount of this cluster generated due to the
	 * annihilation of other clusters.
	 */
	virtual double getGenByAnn();

	/**
	 * This operation retrieves the binding energy for this cluster.
	 * @return An array of the binding energies of this cluster with clusters
	 * of other types as described above.
	 */
	std::vector<double> getBindingEnergy();

	/**
	 * This operation sets the binding energies for this cluster. It expects
	 * the energy vector to be ordered as described above.
	 * @param energies The vector of energies.
	 */
	void setBindingEnergy(std::vector<double> energies);

	/**
	 * This operation retrieves the diffusion factor, D_0, that is used to
	 * calculate the diffusion coefficient for this cluster.
	 * @return The diffusion factor of this cluster
	 */
	double getDiffusionFactor();

	/**
	 * This operation sets the diffusion factor, D_0, that is used to calculate
	 * the diffusion coefficient for this cluster.
	 * @param factor The diffusion factor.
	 */
	void setDiffusionFactor(double factor);

	/**
	 * This operation returns the diffusion coefficient for this cluster and is
	 * calculated from the diffusion factor.
	 * @return The diffusion coefficient.
	 */
	virtual double getDiffusionCoefficient();

};
//end class PSICluster

} /* namespace xolotlCore */
#endif
