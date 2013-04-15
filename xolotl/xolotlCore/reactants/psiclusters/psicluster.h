#ifndef PSICLUSTER_H
#define PSICLUSTER_H

// Includes
#include "reactant.h"

namespace xolotlCore {

// The PSICluster class is a Reactant that is specialized to work for
// simulations of plasma-surface interactions. It provides special routines
// for calculating the total flux due to production and dissociation and
// obtaining the cluster size.
class PSICluster : Reactant {

private:

	// The diffusion factor, D_0, that is used to calculate the diffusion
	// coefficient for this cluster.
	double diffusionFactor;

	// The migration energy for this cluster.
	double migrationEnergy;

	// The binding energy for this cluster.
	double bindingEnergy;

public:

	// Constructor
	PSICluster();

	// Destructor
	~PSICluster();

	// This operation returns the total change in this cluster due to dissociation.
	virtual double getDissociationFlux() = 0;

	// This operation returns the total change in this cluster due to production.
	virtual double getProductionFlux() = 0;

	// This operation returns the total size of the cluster.
	virtual int getSize() = 0;

	// This operation returns the total generation rate of this cluster due to captures.
	virtual double getGenByCapt() = 0;

	// This operation returns the total generation rate of this cluster due to annihilation.
	virtual double getGenByAnn() = 0;

	// This operation retrieves the diffusion factor, D_0, that is used to calculate the diffusion coefficient for this cluster.
	double getDiffusionFactor();

	// This operation retrieves the migration energy for this cluster.
	double getMigrationEnergy();

	// This operation retrieves the binding energy for this cluster.
	double getBindingEnergy();

	// This operation returns the diffusion coefficient for this cluster and is based on the diffusion factor.
	virtual double getDiffusionCoefficient() = 0;

	// This operation sets the diffusion factor, D_0, that is used to calculate the diffusion coefficient for this cluster.
	void setDiffusionFactor(double factor);

	// This operation sets the migration energy for this cluster.
	void setMigrationEnergy(double em);

	// This operation sets the binding energy for this cluster.
	void setBindingEnergy(double eb);

};
//end class PSICluster

} /* namespace xolotlCore */
#endif
