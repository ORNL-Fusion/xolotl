#ifndef DUMMYADVECTIONHANDLER_H
#define DUMMYADVECTIONHANDLER_H

// Includes
#include "SurfaceAdvectionHandler.h"

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile cluster. Here it is a dummy
 * class which means that it is not actually doing anything.
 */
class DummyAdvectionHandler: public SurfaceAdvectionHandler {

public:

	//! The Constructor
	DummyAdvectionHandler() {}

	//! The Destructor
	~DummyAdvectionHandler() {}

	/**
	 * This function initialize the list of clusters that will move through advection. For the
	 * dummy class we don't want any cluster to advect, so this class only clears the vector
	 * and doesn't fill them.
	 *
	 * @param network The network
	 * @param ofill The pointer to the array that will contain the value 1 at the indices
	 * of the advecting clusters
	 */
	void initialize(IReactionNetwork *network, int *ofill) {
		// Clear the index and sink strength vectors
		indexVector.clear();
		sinkStrengthVector.clear();

		// Return now to leave them empty
		return;
	}

};
//end class DummyAdvectionHandler

} /* end namespace xolotlCore */
#endif
