#ifndef FEVCLUSTER_H
#define FEVCLUSTER_H

// Includes
#include "FeCluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of atomic vacancies.
 */
class FeVCluster: public FeCluster {

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
	FeVCluster() = delete;

	/**
	 * The constructor. All FeVClusters must be initialized with a size.
	 *
	 * @param nV the number of atomic vacancies in the cluster
	 * @param registry The performance handler registry
	 */
	FeVCluster(int nV, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			FeCluster(_network, registry, buildName(nV)) {

		// Set the size
		size = nV;
		// Set the typename appropriately
		type = ReactantType::V;

		// Update the composition map
		composition[toCompIdx(Species::V)] = size;

		// Compute the reaction radius
		// It is the same formula for HeV clusters
		reactionRadius = xolotlCore::ironLatticeConstant
				* pow((3.0 * size) / xolotlCore::pi, (1.0 / 3.0)) * 0.5;

		// Bounds on He and V
		heBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		vBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(size),
				static_cast<IReactant::SizeType>(size + 1));

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	FeVCluster(const FeVCluster& other) = delete;

	//! Destructor
	~FeVCluster() {
	}

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation.
	 *
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux() const {
		// Initial declarations
		double flux = FeCluster::getEmissionFlux();

		// Compute the loss to dislocation sinks
		if (size < 5) {
			// k^2 * D * C
			flux += xolotlCore::sinkStrength * diffusionCoefficient
					* concentration;
		}

		return flux;
	}

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getEmissionPartialDerivatives(std::vector<double> & partials) const {
		// Initial declarations
		FeCluster::getEmissionPartialDerivatives(partials);

		// Compute the loss to dislocation sinks
		if (size < 5) {
			// k^2 * D * C
			partials[id - 1] -= xolotlCore::sinkStrength * diffusionCoefficient;
		}

		return;
	}

};
//end class FeVCluster

} /* end namespace xolotlCore */

#endif
