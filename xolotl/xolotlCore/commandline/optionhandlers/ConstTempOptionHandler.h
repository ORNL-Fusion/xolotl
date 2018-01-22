#ifndef CONSTTEMPOPTIONHANDLER_H
#define CONSTTEMPOPTIONHANDLER_H

// Includes
#include <stdlib.h>
#include <TokenizedLineReader.h>
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
					"startTemp <value>                 "
							"The temperature (in Kelvin) will be the constant floating point value specified. "
							"(default = 1000). "
							"\n	                            If two values are given, the second one is interpreted as a gradient."
							"\n	                            (NOTE: Use only ONE temperature option)\n") {
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
	bool handler(IOptions *opt, const std::string& arg) {

		// Set the flag to use constant temperature to true
		opt->setConstTempFlag(true);

		// Build an input stream from the argument string.
		xolotlCore::TokenizedLineReader<std::string> reader;
		auto argSS = std::make_shared < std::istringstream > (arg);
		reader.setInputStream(argSS);

		// Break the argument into tokens.
		auto tokens = reader.loadLine();

		// Set the value for the temperature
		double temp = strtod(tokens[0].c_str(), NULL);
		opt->setConstTemperature(temp);

		// Check if we have another value
		if (tokens.size() > 1) {
			// Set the temperature gradient
			double gradient = strtod(tokens[1].c_str(), NULL);
			opt->setTemperatureGradient(gradient);
		}

		return true;
	}

};
//end class ConstTempOptionHandler

} /* namespace xolotlCore */

#endif
