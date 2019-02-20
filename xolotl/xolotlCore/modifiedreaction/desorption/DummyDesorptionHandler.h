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
			std::vector<double> grid) {
		// Don't do anything
		return;
	}

	/**
	 * This method defines which cluster desorbs at each grid point.
	 *
	 * \see IDesorptionHandler.h
	 */
	void initializeIndex2D(std::vector<int> surfacePos,
			const IReactionNetwork& network, std::vector<double> grid, int ny,
			double hy) {
		// Don't do anything
		return;
	}

	/**
	 * This method defines which cluster desorbs at each grid point.
	 *
	 * \see IDesorptionHandler.h
	 */
	void initializeIndex3D(std::vector<std::vector<int> > surfacePos,
			const IReactionNetwork& network, std::vector<double> grid, int ny,
			double hy, int nz, double hz) {
		// Don't do anything
		return;
	}

};
//end class DummyDesorptionHandler

} /* namespace xolotlCore */
#endif
