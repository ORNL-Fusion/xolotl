#ifndef W100ADVECTIONHANDLER_H
#define W100ADVECTIONHANDLER_H

// Includes
#include "AdvectionHandler.h"

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 */
class W100AdvectionHandler: public AdvectionHandler {

public:

	//! The Constructor
	W100AdvectionHandler() {}

	//! The Destructor
	~W100AdvectionHandler() {}

	/**
	 * This function initialize the list of clusters that will move through advection for a
	 * (100) tungsten material.
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
			if (diffFactor == 0.0) continue;

			// Keep only the helium clusters
			if (cluster->getType() != heType) continue;

			// Get its size
			int heSize = cluster->getSize();

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
			if (sinkStrength == 0.0) continue;

			// Add its index (i) to the vector of indices
			indexVector.push_back(i);

			// Add the sink strength to the vector
			sinkStrengthVector.push_back(sinkStrength);
		}

		return;
	}

};
//end class W100AdvectionHandler

} /* end namespace xolotlCore */
#endif
