#ifndef DENSITYOPTIONHANDLER_H
#define DENSITYOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * DensityOptionHandler handles the depth parameter for bubble bursting.
 */
class DensityOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	DensityOptionHandler() :
			OptionHandler("density",
					"density <rho>             "
							"This option allows the user to set a density in nm-3 "
							"for the number of xenon per volume in a bubble.  \n") {
	}

	/**
	 * The destructor
	 */
	~DensityOptionHandler() {
	}

	/**
	 * This method will set the IOptions density
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The density in nm-3.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Convert to double
		double density = strtod(arg.c_str(), NULL);;
		// Set the number of dimensions
		opt->setDensity(density);

		return true;
	}

};
//end class DensityOptionHandler

} /* namespace xolotlCore */

#endif
