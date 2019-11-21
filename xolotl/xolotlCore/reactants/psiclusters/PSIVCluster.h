#ifndef PSIVCLUSTER_H
#define PSIVCLUSTER_H

// Includes
#include "PSICluster.h"
#include <xolotlPerf.h>
#include <Constants.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of atomic vacancies.
 */
class PSIVCluster: public PSICluster {

private:
	static std::string buildName(IReactant::SizeType nV) {
		std::stringstream nameStream;
		nameStream << "V_" << nV;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSIVCluster() = delete;

	/**
	 * The constructor. All PSIVClusters must be initialized with a size.
	 *
	 * @param nV the number of atomic vacancies in the cluster
	 * @param registry The performance handler registry
	 */
	PSIVCluster(int nV, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSICluster(_network, registry, buildName(nV)) {
		// Set the size
		size = nV;
		// Update the composition map
		composition[toCompIdx(Species::V)] = size;
		// Set the typename appropriately
		type = ReactantType::V;

		// Compute the reaction radius
		// It is the same formula for HeV clusters
		double latticeParam = network.getLatticeParameter();
		reactionRadius = (sqrt(3.0) / 4.0) * latticeParam
				+ pow(
						(3.0 * pow(latticeParam, 3.0) * size)
								/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
				- pow((3.0 * pow(latticeParam, 3.0)) / (8.0 * xolotlCore::pi),
						(1.0 / 3.0));

		// Bounds on He, D, T, and V
		bounds[0] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		bounds[1] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		bounds[2] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		bounds[3] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(size),
				static_cast<IReactant::SizeType>(size + 1));

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	PSIVCluster(const PSIVCluster& other) = delete;

	//! Destructor
	~PSIVCluster() {
	}

	/**
	 * Add grid points to the vector of diffusion coefficients or remove
	 * them if the value is negative.
	 *
	 * @param i The number of grid point to add or remove
	 */
	void addGridPoints(int i) override {
		if (diffusionFactor > 0.0) {
			Reactant::addGridPoints(i);
		}

		// Don't do anything
		return;
	}

	/**
	 * This operation sets the temperature at which the reactant currently
	 * exists. Temperature-dependent quantities are recomputed when this
	 * operation is called, so the temperature should always be set first.
	 *
	 * @param temp The new cluster temperature
	 * @param i The location on the grid
	 */
	void setTemperature(double temp, int i) override {
		if (diffusionFactor > 0.0) {
			Reactant::setTemperature(temp, i);
		}

		// Don't do anything
		return;
	}

	/**
	 * This operation returns the diffusion coefficient for this reactant and is
	 * calculated from the diffusion factor.
	 *
	 * @param i The position on the grid
	 * @return The diffusion coefficient
	 */
	double getDiffusionCoefficient(int i) const override {
		if (diffusionFactor > 0.0) {
			return Reactant::getDiffusionCoefficient(i);
		}

		return 0.0;
	}

};
//end class PSIVCluster

} /* end namespace xolotlCore */

#endif
