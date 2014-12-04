#ifndef CONSTTEMPOPTIONHANDLER_H
#define CONSTTEMPOPTIONHANDLER_H

// Includes
#include <stdlib.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * ConstTempOptionHandler handles the case where we want to use a constant temperature.
 */
class ConstTempOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	ConstTempOptionHandler() :
			OptionHandler("startTemp",
					"startTemp <value>           "
							"The temperature (in Kelvin) will be the constant floating point value specified. "
							"(default = 1000)"
							"\n	                      (NOTE: Use only ONE temperature option)") {
	}

	/**
	 * The destructor
	 */
	~ConstTempOptionHandler() {
	}

	/**
	 * This method will set the IOptions constTempFlag and constTemperature
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The value for the temperature.
	 */
	bool handler(IOptions *opt, std::string arg) {

		// Set the flag to use constant temperature to true
		opt->setConstTempFlag(true);

		// Set the value for the temperature
		double temp = strtod(arg.c_str(), NULL);
		opt->setConstTemperature(temp);

		return true;
	}

};
//end class ConstTempOptionHandler

} /* namespace xolotlCore */

#endif
