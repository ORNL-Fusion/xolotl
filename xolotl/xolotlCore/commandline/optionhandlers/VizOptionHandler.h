#ifndef VIZOPTIONHANDLER_H
#define VIZOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * VizOptionHandler handles the choice of handlers for the visualization infrastructure.
 */
class VizOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	VizOptionHandler() :
		OptionHandler("vizHandler",
				"vizHandler {std,dummy}            "
				"Which set of handlers to use for the visualization. (default = dummy)\n") {}

	/**
	 * The destructor
	 */
	~VizOptionHandler() {
	}

	/**
	 * This method will set the IOptions vizStandardHandlersFlag
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The argument for the flag.
	 */
	bool handler(IOptions *opt, std::string arg) {
		// Determine the type of handlers we are being asked to use
		if (arg == "std") {
			opt->setVizStandardHandlers(true);
		}
		else if (arg == "dummy") {
			opt->setVizStandardHandlers(false);
		}
		else {
			std::cerr << "Options: unrecognized argument in the visualization option handler: " << arg << std::endl;
			opt->showHelp(std::cerr);
			opt->setShouldRunFlag(false);
			opt->setExitCode(EXIT_FAILURE);
			return false;
		}

		return true;
	}

};
//end class VizOptionHandler

} /* namespace xolotlCore */

#endif
