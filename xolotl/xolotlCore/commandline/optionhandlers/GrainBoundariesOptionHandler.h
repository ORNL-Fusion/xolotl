#ifndef GRAINBOUNDARIESOPTIONHANDLER_H
#define GRAINBOUNDARIESOPTIONHANDLER_H

// Includes
#include <string.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * GrainBoundariesOptionHandler handles the option to generate
 * grain boundaries in the material.
 */
class GrainBoundariesOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	GrainBoundariesOptionHandler() :
		OptionHandler("grain",
				"grain                         "
				"This option allows the user to add GB in the X, Y, or Z directions.\n"
				"                                    To do so, simply write the direction followed "
				"by the distance in nm, for instance: X 3.0 Z 2.5 Z 10.0 .\n") {}

	/**
	 * The destructor
	 */
	~GrainBoundariesOptionHandler() {
	}

	/**
	 * This method will set the IOptions gbString
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The list of PETSc options.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
    	// Set the material name
    	opt->setGbString(arg);
    	return true;
	}

};
//end class GrainBoundariesOptionHandler

} /* namespace xolotlCore */

#endif
