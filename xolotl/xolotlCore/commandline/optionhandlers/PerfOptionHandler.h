#ifndef PERFOPTIONHANDLER_H
#define PERFOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * PerfOptionHandler handles the choice of handlers for the performance infrastructure.
 */
class PerfOptionHandler: public OptionHandler {
protected:

	/**
	 * The default constructor
	 */
	PerfOptionHandler() :
			OptionHandler() {
	}
	;

public:

	/**
	 * The constructor to use.
	 * @param keyName The name for the key.
	 * @param msg The help message.
	 */
	PerfOptionHandler(std::string keyName, std::string msg) :
			OptionHandler(keyName, msg) {
	}

	/**
	 * The destructor
	 */
	~PerfOptionHandler() {
	}

	/**
	 * This method will set the IOptions perfStandardHandlersFlag
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The argument for the flag.
	 */
	bool handler(IOptions *opt, std::string arg) {
		// Determine the type of handlers we are being asked to use
		if (arg == "std") {
			opt->setPerfStandardHandlers(true);
		}
		else if (arg == "dummy") {
			opt->setPerfStandardHandlers(false);
		}
		else {
			std::cerr << "Options: unrecognized argument in the performance option handler: " << arg << std::endl;
			opt->showHelp(std::cerr);
			opt->setShouldRunFlag(false);
			opt->setExitCode(EXIT_FAILURE);
			return false;
		}

		return true;
	}

};
//end class PerfOptionHandler

} /* namespace xolotlCore */

#endif
