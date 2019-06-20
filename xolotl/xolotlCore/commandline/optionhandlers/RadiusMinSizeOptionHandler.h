#ifndef RADIUSMINSIZEOPTIONHANDLER_H
#define RADIUSMINSIZEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * RadiusMinSizeOptionHandler handles the number of dimensions option.
 */
class RadiusMinSizeOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	RadiusMinSizeOptionHandler() :
			OptionHandler("radiusSize",
					"radiusSize <minSize>              "
							"This option allows the user a minimum size for the computation for the average radius (default is 0).  \n") {
	}

	/**
	 * The destructor
	 */
	~RadiusMinSizeOptionHandler() {
	}

	/**
	 * This method will set the IOptions resoMinSize
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The minimum size.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Convert to integer
		int size = strtol(arg.c_str(), NULL, 10);
		// Set the number of dimensions
		opt->setRadiusMinSize(size);

		return true;
	}

};
//end class RadiusMinSizeOptionHandler

} /* namespace xolotlCore */

#endif
