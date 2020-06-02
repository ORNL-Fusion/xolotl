#ifndef W211TRAPMUTATIONHANDLER_H
#define W211TRAPMUTATIONHANDLER_H

#include <xolotl/core/modified/TrapMutationHandler.h>

namespace xolotl {
namespace core {
namespace modified {

/**
 * This class realizes the ITrapMutationHandler interface responsible for the modified
 * trap-mutation of small helium clusters close to the surface for a (211) oriented
 * tungsten material.
 */
class W211TrapMutationHandler: public TrapMutationHandler {
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
			depthVec = {0.5, 0.8, 1.0, 1.0, 1.0, 1.3, 1.5};
			sizeVec = {1, 1, 1, 2, 2, 2, 2};

			// He1 desorpts with 64%
			desorp = Desorption(1, 0.64);
		}
		else {
			depthVec = {0.5, 0.8, 1.0, 1.0, 1.3, 1.3, 1.2};
			sizeVec = {1, 1, 1, 2, 1, 2, 3};

			// He1 desorpts with 59%
			desorp = Desorption(1, 0.59);
		}

		return;
	}

public:

	/**
	 * The constructor
	 */
	W211TrapMutationHandler() {}

	/**
	 * The Destructor
	 */
	~W211TrapMutationHandler() {}

};
//end class W211TrapMutationHandler

} /* namespace modified */
} /* namespace core */
} /* namespace xolotl */

#endif
