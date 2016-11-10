#ifndef VOIDPORTIONOPTIONHANDLER_H
#define VOIDPORTIONOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * VoidPortionOptionHandler handles the use of the void portion at the start of the
 * simulation option. It represents the room that is left for the surface to grow.
 */
class VoidPortionOptionHandler : public OptionHandler {
public:

	/**
	 * The default constructor
	 */
    VoidPortionOptionHandler() :
    	OptionHandler("voidPortion",
    			"voidPortion <value>                 "
    			"The value (in %) of the void portion at the start of the simulation.") {}

	/**
	 * The destructor
	 */
    ~VoidPortionOptionHandler() {}

    /**
     * This method will set the IOptions voidPortion to the value given
     * as the argument.
     *
     * @param opt The pointer to the option that will be modified.
     * @param arg The value for the void portion.
     */
    bool handler(IOptions *opt, const std::string& arg) {
    	// Set the value
    	double portion = strtod(arg.c_str(), NULL);

    	// Check the given value makes sense
    	if (portion > 100.0 || portion < 0.0) {
    		std::cerr << "Options: wrong value for the void portion option handler: " << arg << std::endl;
    		opt->showHelp(std::cerr);
    		opt->setShouldRunFlag(false);
    		opt->setExitCode(EXIT_FAILURE);
    		return false;
    	}

    	opt->setVoidPortion(portion);
    	return true;
    }

};//end class VoidPortionOptionHandler

} /* namespace xolotlCore */

#endif
