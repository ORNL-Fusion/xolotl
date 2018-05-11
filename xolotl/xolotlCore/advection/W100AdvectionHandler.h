#ifndef W100ADVECTIONHANDLER_H
#define W100ADVECTIONHANDLER_H

// Includes
#include "SurfaceAdvectionHandler.h"
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 */
class W100AdvectionHandler: public SurfaceAdvectionHandler {

public:

	//! The Constructor
	W100AdvectionHandler() :
			SurfaceAdvectionHandler() {
	}

	//! The Destructor
	~W100AdvectionHandler() {
	}

	/**
	 * This function initialize the list of clusters that will move through advection for a
	 * (100) tungsten material.
	 *
	 * @param network The network
	 * @param ofillMap Map of connectivity for advecting clusters.
	 * of the advecting clusters
	 */
	// TODO this is nearly identical for the W100, W110, W111, and W211
	// cases.  Factor the identical parts to a base class, and only
	// have these classes differ in the sinkStrength identification.
	void initialize(const IReactionNetwork& network,
                    IReactionNetwork::SparseFillMap& ofillMap) override {
		// Get all the reactants and their number
		int dof = network.getDOF();

		// Clear the index and sink strength vectors
		advectingClusters.clear();
		sinkStrengthVector.clear();

		// Consider each reactant.
		for (IReactant const& currReactant : network.getAll()) {

			auto const& cluster = static_cast<PSICluster const&>(currReactant);

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
				sinkStrength = 2.28e-3;
				break;
			case 2:
				sinkStrength = 5.06e-3;
				break;
			case 3:
				sinkStrength = 7.26e-3;
				break;
			case 4:
				sinkStrength = 15.87e-3;
				break;
			case 5:
				sinkStrength = 16.95e-3;
				break;
			case 6:
				sinkStrength = 27.16e-3;
				break;
			case 7:
				sinkStrength = 35.56e-3;
				break;
			}

			// If the sink strength is still 0.0, this cluster is not advecting
			if (xolotlCore::equal(sinkStrength, 0.0))
				continue;

			// Add it to our collection of advecting clusters.
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
//end class W100AdvectionHandler

} /* end namespace xolotlCore */
#endif
