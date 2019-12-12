#ifndef DUMMYDESORPTIONHANDLER_H
#define DUMMYDESORPTIONHANDLER_H

// Includes
#include <DesorptionHandler.h>

namespace xolotlCore {

/**
 * This class realizes IDesorptionHandler interface, responsible for the desorption
 *  of hydrogen clusters at the surface. Here it is a dummy class,
 * meaning that it should not do anything.
 */
class DummyDesorptionHandler: public DesorptionHandler {

public:

	/**
	 * The constructor
	 */
	DummyDesorptionHandler() {
		nCluster = 0;
		kRecombination = 0.0;
		equilibriumConc = 0.0;
	}

	/**
	 * The destructor
	 */
	~DummyDesorptionHandler() {
	}

	/**
	 * This method defines which cluster desorbs at each grid point.
	 *
	 * \see IDesorptionHandler.h
	 */
	void initializeIndex1D(int surfacePos, const IReactionNetwork& network,
			int nx, int xs) {
		// Don't do anything
		return;
	}

	/**
	 * This method defines which cluster desorbs at each grid point.
	 *
	 * \see IDesorptionHandler.h
	 */
	void initializeIndex2D(std::vector<int> surfacePos,
			const IReactionNetwork& network, int nx, int xs, int ny = 1,
			int ys = 0) {
		// Don't do anything
		return;
	}

	/**
	 * This method defines which cluster desorbs at each grid point.
	 *
	 * \see IDesorptionHandler.h
	 */
	void initializeIndex3D(std::vector<std::vector<int> > surfacePos,
			const IReactionNetwork& network, int nx, int xs, int ny = 1,
			int ys = 0, int nz = 1, int zs = 0) {
		// Don't do anything
		return;
	}

};
//end class DummyDesorptionHandler

} /* namespace xolotlCore */
#endif
