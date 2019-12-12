#ifndef W111ADVECTIONHANDLER_H
#define W111ADVECTIONHANDLER_H

// Includes
#include "SurfaceAdvectionHandler.h"
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 */
class W111AdvectionHandler: public SurfaceAdvectionHandler {

public:

	//! The Constructor
	W111AdvectionHandler() :
			SurfaceAdvectionHandler() {
	}

	//! The Destructor
	~W111AdvectionHandler() {
	}

	/**
	 * This function initialize the list of clusters that will move through advection for a
	 * (111) tungsten material.
	 *
	 * @param network The network
	 * @param ofillMap Map of connectivity for advecting clusters.
	 * of the advecting clusters
	 */
	void initialize(const IReactionNetwork& network,
                    IReactionNetwork::SparseFillMap& ofillMap) override {

		int dof = network.getDOF();

		// Clear the index and sink strength vectors
		advectingClusters.clear();
		sinkStrengthVector.clear();

		// Consider each reactant.
		for (IReactant const& currReactant : network.getAll()) {

			auto const& cluster = static_cast<IReactant const&>(currReactant);

			// Get its diffusion coefficient
			double diffFactor = cluster.getDiffusionFactor();

			// Don't do anything if the diffusion factor is 0.0
			if (xolotlCore::equal(diffFactor, 0.0))
				continue;

			// Keep only the helium clusters
			if (cluster.getType() != ReactantType::He)
				continue;

			// Get its size
			int heSize = cluster.getSize();

			// Switch on the size to get the sink strength (in eV.nm3)
			double sinkStrength = 0.0;
			switch (heSize) {
			case 1:
				sinkStrength = 3.65e-3;
				break;
			case 2:
				sinkStrength = 6.40e-3;
				break;
			case 3:
				sinkStrength = 16.38e-3;
				break;
			case 4:
				sinkStrength = 9.84e-3;
				break;
			case 5:
				sinkStrength = 44.40e-3;
				break;
			case 6:
				sinkStrength = 52.12e-3;
				break;
			case 7:
				sinkStrength = 81.57e-3;
				break;
			}

			// If the sink strength is still 0.0, this cluster is not advecting
			if (xolotlCore::equal(sinkStrength, 0.0))
				continue;

			// Add it to the collection of advecting clusters.
			advectingClusters.emplace_back(currReactant);

			// Add the sink strength to the vector
			sinkStrengthVector.push_back(sinkStrength);

			// Set the off-diagonal part for the Jacobian to 1
			// Get its id
			int index = cluster.getId() - 1;
			// Set the ofill value to 1 for this cluster
            ofillMap[index].emplace_back(index);
		}

		return;
	}

};
//end class W111AdvectionHandler

} /* end namespace xolotlCore */
#endif
