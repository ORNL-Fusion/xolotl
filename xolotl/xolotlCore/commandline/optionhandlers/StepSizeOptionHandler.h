#ifndef STEPSIZEOPTIONHANDLER_H
#define STEPSIZEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * StepSizeOptionHandler handles the choice of handlers for the visualization infrastructure.
 */
class StepSizeOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	StepSizeOptionHandler() :
		OptionHandler("stepSize",
				"stepSize <value>      "
				"The size of the step for the spatial grid in nm. (default = 1.0nm)") {}

	/**
	 * The destructor
	 */
	~StepSizeOptionHandler() {
	}

	/**
	 * This method will set the IOptions stepSize
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The argument for the flag.
	 */
	bool handler(IOptions *opt, std::string arg) {
    	// Set the value for the step
    	double step = strtod(arg.c_str(), NULL);
    	opt->setStepSize(step);

		return true;
	}

};
//end class StepSizeOptionHandler

} /* namespace xolotlCore */

#endif
