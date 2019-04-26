#ifndef BURSTINGMINSIZEOPTIONHANDLER_H
#define BURSTINGMINSIZEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * BurstingMinSizeOptionHandler handles the minimum size for bubble bursting.
 */
class BurstingMinSizeOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	BurstingMinSizeOptionHandler() :
			OptionHandler("burstingMin", "burstingMin <size>                "
					"This option allows the user to set a minimum size "
					"for the bubble bursting.  \n") {
	}

	/**
	 * The destructor
	 */
	~BurstingMinSizeOptionHandler() {
	}

	/**
	 * This method will set the IOptions burstingDepth
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The depth in nm.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Convert to double
		int size = strtol(arg.c_str(), NULL, 10);
		// Set the number of dimensions
		opt->setBurstingSize(size);

		return true;
	}

};
//end class BurstingMinSizeOptionHandler

} /* namespace xolotlCore */

#endif
