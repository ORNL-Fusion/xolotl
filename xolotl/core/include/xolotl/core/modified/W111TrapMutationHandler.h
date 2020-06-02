#ifndef W111TRAPMUTATIONHANDLER_H
#define W111TRAPMUTATIONHANDLER_H

#include <xolotl/core/modified/TrapMutationHandler.h>

namespace xolotl {
namespace core {
namespace modified {

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
	 *
	 * @param temp The temperature of the system
	 */
	void initializeDepthSize(double temp) {
		// Switch values depending on the temperature
		if (temp < 1066.5) {
			depthVec = {0.6, 0.8, 1.1, 1.1, 1.2, 1.3, 1.3};
			sizeVec = {1, 1, 1, 1, 1, 1, 2};

			// He1 desorpts with 61%
			desorp = Desorption(1, 0.61);
		}
		else {
			depthVec = {0.6, 0.8, 1.1, 1.1, 1.1, 1.1, 1.1};
			sizeVec = {1, 1, 1, 1, 1, 1, 2};

			// He1 desorpts with 35%
			desorp = Desorption(1, 0.35);
		}

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

} /* namespace modified */
} /* namespace core */
} /* namespace xolotl */

#endif
