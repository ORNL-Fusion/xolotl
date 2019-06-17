#ifndef RESOMINSIZEOPTIONHANDLER_H
#define RESOMINSIZEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * ResoMinSizeOptionHandler handles the number of dimensions option.
 */
class ResoMinSizeOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	ResoMinSizeOptionHandler() :
			OptionHandler("resoSize",
					"resoSize <minSize>                "
							"This option allows the user a minimum size for the re-solution (default is 0).  \n") {
	}

	/**
	 * The destructor
	 */
	~ResoMinSizeOptionHandler() {
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
		opt->setResoMinSize(size);

		return true;
	}

};
//end class ResoMinSizeOptionHandler

} /* namespace xolotlCore */

#endif
