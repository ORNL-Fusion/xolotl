#ifndef W110TRAPMUTATIONHANDLER_H
#define W110TRAPMUTATIONHANDLER_H

#include <xolotl/core/modified/TrapMutationHandler.h>

namespace xolotl
{
namespace core
{
namespace modified
{
/**
 * This class realizes the ITrapMutationHandler interface responsible for the
 * modified trap-mutation of small helium clusters close to the surface for a
 * (110) oriented tungsten material.
 */
class W110TrapMutationHandler : public TrapMutationHandler
{
private:
	/**
	 * \see TrapMutationHandler.h
	 */
	void
	initializeDepthSize(double temp)
	{
		// Switch values depending on the temperature
		if (temp < 1066.5) {
			depthVec = {-0.1, 0.7, 0.9, 0.9, 0.9, 1.1, 1.1};
			sizeVec = {0, 1, 1, 1, 1, 2, 2};

			// He2 desorpts with 31%
			desorp = Desorption(2, 0.31);
		}
		else {
			depthVec = {-0.1, 0.7, 0.9, 0.9, 0.9, 1.1, 1.1};
			sizeVec = {0, 1, 1, 1, 1, 1, 2};

			// He2 desorpts with 32%
			desorp = Desorption(2, 0.32);
		}

		return;
	}

public:
	/**
	 * The constructor
	 */
	W110TrapMutationHandler()
	{
	}

	/**
	 * The Destructor
	 */
	~W110TrapMutationHandler()
	{
	}
};
// end class W110TrapMutationHandler

} /* namespace modified */
} /* namespace core */
} /* namespace xolotl */

#endif
