#ifndef W211ADVECTIONHANDLER_H
#define W211ADVECTIONHANDLER_H

// Includes
#include "SurfaceAdvectionHandler.h"
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 */
class W211AdvectionHandler: public SurfaceAdvectionHandler {

public:

	//! The Constructor
	W211AdvectionHandler() :
			SurfaceAdvectionHandler() {
	}

	//! The Destructor
	~W211AdvectionHandler() {
	}

	/**
	 * This function initialize the list of clusters that will move through advection for a
	 * (211) tungsten material.
	 *
	 * @param network The network
	 * @param ofillMap Map of connectivity for advecting clusters.
	 * of the advecting clusters
	 */
	void initialize(experimental::IReactionNetwork& network,
			experimental::IReactionNetwork::SparseFillMap& ofillMap) override {
		// Clear the index and sink strength vectors
		advectingClusters.clear();
		sinkStrengthVector.clear();

		using NetworkType =
		experimental::PSIReactionNetwork<experimental::PSIFullSpeciesList>;
		auto psiNetwork = dynamic_cast<NetworkType*>(&network);

		// Initialize the composition
		NetworkType::Composition comp;
		for (auto i : psiNetwork->getSpeciesRange()) {
			comp[i] = 0;
		}

		// Loop on helium clusters from size 1 to 7
		for (std::size_t i = 0; i < 7; i++) {
			comp[NetworkType::Species::He] = i;
			auto cluster = psiNetwork->findCluster(comp, plsm::onHost);

			// Get its diffusion coefficient
			double diffFactor = cluster.getDiffusionFactor();

			// Don't do anything if the diffusion factor is 0.0
			if (xolotlCore::equal(diffFactor, 0.0))
				continue;

			// Switch on the size to get the sink strength (in eV.nm3)
			double sinkStrength = 0.0;
			switch (i) {
			case 1:
				sinkStrength = 1.49e-3;
				break;
			case 2:
				sinkStrength = 3.69e-3;
				break;
			case 3:
				sinkStrength = 12.34e-3;
				break;
			case 4:
				sinkStrength = 14.11e-3;
				break;
			case 5:
				sinkStrength = 19.14e-3;
				break;
			case 6:
				sinkStrength = 35.77e-3;
				break;
			case 7:
				sinkStrength = 67.65e-3;
				break;
			}

			// If the sink strength is still 0.0, this cluster is not advecting
			if (xolotlCore::equal(sinkStrength, 0.0))
				continue;

			// Get its id
			auto index = cluster.getId();
			// Add it to our collection of advecting clusters.
			advectingClusters.emplace_back(index);

			// Add the sink strength to the vector
			sinkStrengthVector.push_back(sinkStrength);

			// Set the off-diagonal part for the Jacobian to 1
			// Set the ofill value to 1 for this cluster
			ofillMap[index].emplace_back(index);
		}

		return;
	}

};
//end class W211AdvectionHandler

} /* end namespace xolotlCore */
#endif
