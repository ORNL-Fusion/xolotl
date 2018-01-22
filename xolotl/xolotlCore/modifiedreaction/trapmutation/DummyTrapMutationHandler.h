#ifndef DUMMYTRAPMUTATIONHANDLER_H
#define DUMMYTRAPMUTATIONHANDLER_H

#include "TrapMutationHandler.h"

namespace xolotlCore {

/**
 * This class realizes the ITrapMutationHandler interface responsible for the modified
 * trap-mutation of small helium clusters close to the surface. Here it is a dummy class,
 * meaning that it should not do anything.
 */
class DummyTrapMutationHandler: public TrapMutationHandler {
private:

	/**
	 * Method initializing the depth vector, the size vector,
	 * and desorption information.
	 *
	 * @param temp The temperature of the system
	 */
	void initializeDepthSize(double temp) {
		// Initialize the vectors
		depthVec = {-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1};
		sizeVec = {0, 0, 0, 0, 0, 0, 0};

		// And don't do anything else
		return;
	}

public:

	/**
	 * The constructor
	 */
	DummyTrapMutationHandler() {}

	/**
	 * The Destructor
	 */
	~DummyTrapMutationHandler() {}

};
//end class DummyTrapMutationHandler

}

#endif
