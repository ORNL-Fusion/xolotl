#ifndef W110ADVECTIONHANDLER_H
#define W110ADVECTIONHANDLER_H

// Includes
#include "AdvectionHandler.h"

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile cluster.
 */
class W110AdvectionHandler: public AdvectionHandler {

public:

	//! The Constructor
	W110AdvectionHandler() {}

	//! The Destructor
	~W110AdvectionHandler() {}

	/**
	 * The off-diagonal part of the Jacobian is already initialized by the diffusion handler.
	 * This function initialize the list of clusters that will move through advection for a
	 * (110) tungsten material.
	 *
	 * @param network The network
	 */
	void initialize(PSIClusterReactionNetwork *network) {
		// Get all the reactant
		auto reactants = network->getAll();
		int size = reactants->size();

		// Clear the index and sink strength vectors
		indexVector.clear();
		sinkStrengthVector.clear();

		// Loop on them
		for (int i = 0; i < size; i++) {
			// Get the i-th cluster
			auto cluster = (PSICluster *) reactants->at(i);
			// Get its diffusion coefficient
			double diffFactor = cluster->getDiffusionFactor();

			// Don't do anything if the diffusion factor is 0.0
			if (diffFactor == 0.0) continue;

			// Only helium clusters
			if (cluster->getType() != heType) continue;

			// Get its size
			int size = cluster->getSize();

			// Switch on it to get the sink strength (in eV.nm3)
			double sinkStrength = 0.0;
			switch (size) {
				case 1:
					sinkStrength = 0.92e-3;
					break;
				case 2:
					sinkStrength = 1.48e-3;
					break;
				case 3:
					sinkStrength = 6.73e-3;
					break;
				case 4:
					sinkStrength = 6.18e-3;
					break;
				case 5:
					sinkStrength = 33.61e-3;
					break;
				case 6:
					sinkStrength = 37.58e-3;
					break;
				case 7:
					sinkStrength = 41.90e-3;
					break;
			}

			// If the sink strength is still 0.0, this cluster is not advecting
			if (sinkStrength == 0.0) continue;

			// Add it's index (i) to the vector of indices
			indexVector.push_back(i);

			// Add the sink strength to the vector
			sinkStrengthVector.push_back(sinkStrength);
		}

		return;
	}

};
//end class W110AdvectionHandler

} /* end namespace xolotlCore */
#endif
