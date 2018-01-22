#ifndef W100TRAPMUTATIONHANDLER_H
#define W100TRAPMUTATIONHANDLER_H

#include "TrapMutationHandler.h"

namespace xolotlCore {

/**
 * This class realizes the ITrapMutationHandler interface responsible for the modified
 * trap-mutation of small helium clusters close to the surface for a (100) oriented
 * tungsten material.
 */
class W100TrapMutationHandler: public TrapMutationHandler {
private:

	/**
	 * Method initializing the depth vector, the size vector,
	 * and desorption information.
	 *
	 * @param temp The temperature of the system
	 */
	void initializeDepthSize(double temp) {
		// Switch values depending on the temperature
		if (temp < 1066.5) {
			depthVec = {-0.1, 0.5, 0.6, 0.6, 0.8, 0.8, 0.8};
			sizeVec = {0, 1, 1, 1, 1, 2, 2};

			// He2 desorpts with 4%
			desorp = Desorption(2, 0.04);
		}
		else {
			depthVec = {-0.1, 0.5, 0.6, 0.8, 0.6, 0.8, 0.8};
			sizeVec = {0, 1, 1, 1, 2, 2, 2};

			// He2 desorpts with 19%
			desorp = Desorption(2, 0.19);
		}

		return;
	}

public:

	/**
	 * The constructor
	 */
	W100TrapMutationHandler() {}

	/**
	 * The Destructor
	 */
	~W100TrapMutationHandler() {}

};
//end class W100TrapMutationHandler

}

#endif
