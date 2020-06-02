#ifndef SIGMA3TRAPMUTATIONHANDLER_H
#define SIGMA3TRAPMUTATIONHANDLER_H

namespace xolotl {
namespace core {
namespace modified {

/**
 * This class implements the modified trap mutation occurring near a sigma 3
 * grain boundary. It helps the main trap mutation handler filling its index vector.
 */
class Sigma3TrapMutationHandler {
private:

	/** The vector containing the different distances for the modified trap-mutation
	 * associated with the GB
	 */
	std::vector<double> distanceVec = {0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5};

	//! The vector containing the different vacancy size for the modified trap-mutation
	std::vector<int> sizeVec = {0, 0, 0, 1, 1, 1, 1};

public:

	/**
	 * The constructor
	 */
	Sigma3TrapMutationHandler() {}

	/**
	 * The Destructor
	 */
	~Sigma3TrapMutationHandler() {}

	/**
	 * Returns the distance vector.
	 */
	std::vector<double> getDistanceVector() {return distanceVec;}

	/**
	 * Returns the size vector.
	 */
	std::vector<int> getSizeVector() {return sizeVec;}

};
//end class Sigma3TrapMutationHandler

} /* namespace modified */
} /* namespace core */
} /* namespace xolotl */

#endif
