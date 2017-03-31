#ifndef DIMENSIONSOPTIONHANDLER_H
#define DIMENSIONSOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * DimensionsOptionHandler handles the number of dimensions option.
 */
class DimensionsOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	DimensionsOptionHandler() :
			OptionHandler("dimensions",
					"dimensions <dimNumber>            "
							"This option allows the user to switch between a 0, 1, 2, or 3 "
							"dimensions simulation.  \n") {
	}

	/**
	 * The destructor
	 */
	~DimensionsOptionHandler() {
	}

	/**
	 * This method will set the IOptions dimensionNumber
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The number of dimensions.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Convert to integer
		int number = strtol(arg.c_str(), NULL, 10);
		// Set the number of dimensions
		opt->setDimensionNumber(number);

		return true;
	}

};
//end class DimensionsOptionHandler

} /* namespace xolotlCore */

#endif
