#ifndef W110TRAPMUTATIONHANDLER_H
#define W110TRAPMUTATIONHANDLER_H

#include "TrapMutationHandler.h"

namespace xolotlCore {

/**
 * This class realizes the ITrapMutationHandler interface responsible for the modified
 * trap-mutation of small helium clusters close to the surface for a (110) oriented
 * tungsten material.
 */
class W110TrapMutationHandler: public TrapMutationHandler {
private:

	/**
	 * Method initializing the depth vector, the size vector,
	 * and desorption information.
	 */
	void initializeDepthSize() {
		depthVec = {-0.1, 0.7, 0.9, 0.9, 0.9, 1.1, 1.1};
		sizeVec = {0, 1, 1, 1, 1, 1, 2};

		// He2 desorpts with 32%
		desorp = Desorption(2, 0.32);

		return;
	}

public:

	/**
	 * The constructor
	 */
	W110TrapMutationHandler() {}

	/**
	 * The Destructor
	 */
	~W110TrapMutationHandler() {}

};
//end class W110TrapMutationHandler

}

#endif
