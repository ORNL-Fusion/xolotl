#ifndef BURSTINGDEPTHOPTIONHANDLER_H
#define BURSTINGDEPTHOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * BurstingDepthOptionHandler handles the depth parameter for bubble bursting.
 */
class BurstingDepthOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	BurstingDepthOptionHandler() :
			OptionHandler("burstingDepth",
					"burstingDepth <depth>             "
							"This option allows the user to set a depth in nm "
							"for the bubble bursting.  \n") {
	}

	/**
	 * The destructor
	 */
	~BurstingDepthOptionHandler() {
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
		double depth = strtod(arg.c_str(), NULL);;
		// Set the number of dimensions
		opt->setBurstingDepth(depth);

		return true;
	}

};
//end class BurstingDepthOptionHandler

} /* namespace xolotlCore */

#endif
