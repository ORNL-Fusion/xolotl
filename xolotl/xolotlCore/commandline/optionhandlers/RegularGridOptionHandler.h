#ifndef REGULARGRIDOPTIONHANDLER_H
#define REGULARGRIDOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * RegularGridOptionHandler handles the choice for the grid in the x direction.
 * It will be either regular or not.
 */
class RegularGridOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	RegularGridOptionHandler() :
			OptionHandler("regularGrid", "regularGrid {yes,  no}            "
					"Will the grid be regularly spaced in the x direction?\n") {
	}

	/**
	 * The destructor
	 */
	~RegularGridOptionHandler() {
	}

	/**
	 * This method will set the IOptions useRegularGridFlag
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The argument for the flag.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Determine the type of handlers we are being asked to use
		if (arg == "yes") {
			opt->setRegularXGrid(true);
		} else if (arg == "no") {
			opt->setRegularXGrid(false);
		} else if (arg == "cheby") {
			opt->setChebyshevGrid(true);
		} else {
			std::cerr
					<< "Options: unrecognized argument in the regular grid option handler: "
					<< arg << std::endl;
			opt->showHelp(std::cerr);
			opt->setShouldRunFlag(false);
			opt->setExitCode(EXIT_FAILURE);
			return false;
		}

		return true;
	}

};
//end class RegularGridOptionHandler

} /* namespace xolotlCore */

#endif
