#ifndef W111ADVECTIONHANDLER_H
#define W111ADVECTIONHANDLER_H

// Includes
#include "AdvectionHandler.h"
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 */
class W111AdvectionHandler: public AdvectionHandler {

public:

	//! The Constructor
	W111AdvectionHandler() {}

	//! The Destructor
	~W111AdvectionHandler() {}

	/**
	 * This function initialize the list of clusters that will move through advection for a
	 * (111) tungsten material.
	 *
	 * @param network The network
	 */
	void initialize(PSIClusterReactionNetwork *network) {
		// Get all the reactants and their number
		auto reactants = network->getAll();
		int size = reactants->size();

		// Clear the index and sink strength vectors
		indexVector.clear();
		sinkStrengthVector.clear();

		// Loop on all the reactants
		for (int i = 0; i < size; i++) {
			// Get the i-th cluster
			auto cluster = (PSICluster *) reactants->at(i);
			// Get its diffusion coefficient
			double diffFactor = cluster->getDiffusionFactor();

			// Don't do anything if the diffusion factor is 0.0
			if (xolotlCore::equal(diffFactor, 0.0)) continue;

			// Keep only the helium clusters
			if (cluster->getType() != heType) continue;

			// Get its size
			int heSize = cluster->getSize();

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
			if (xolotlCore::equal(sinkStrength, 0.0)) continue;

			// Add its index (i) to the vector of indices
			indexVector.push_back(i);

			// Add the sink strength to the vector
			sinkStrengthVector.push_back(sinkStrength);
		}

		return;
	}

};
//end class W111AdvectionHandler

} /* end namespace xolotlCore */
#endif
