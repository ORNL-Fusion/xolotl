#ifndef W111TRAPMUTATIONHANDLER_H
#define W111TRAPMUTATIONHANDLER_H

#include "TrapMutationHandler.h"

namespace xolotlCore {

/**
 * This class realizes the ITrapMutationHandler interface responsible for the modified
 * trap-mutation of small helium clusters close to the surface for a (111) oriented
 * tungsten material.
 */
class W111TrapMutationHandler: public TrapMutationHandler {
private:

	/**
	 * Method initializing the depth vector, the size vector,
	 * and desorption information.
	 */
	void initializeDepthSize() {
		depthVec = {0.6, 0.8, 1.1, 1.1, 1.1, 1.1, 1.1};
		sizeVec = {1, 1, 1, 1, 1, 1, 2};

		// He1 desorpts with 35%
		desorp = Desorption(1, 0.35);

		return;
	}

public:

	/**
	 * The constructor
	 */
	W111TrapMutationHandler() {}

	/**
	 * The Destructor
	 */
	~W111TrapMutationHandler() {}

};
//end class W111TrapMutationHandler

}

#endif
