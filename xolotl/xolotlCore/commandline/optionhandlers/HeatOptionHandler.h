#ifndef HEATOPTIONHANDLER_H
#define HEATOPTIONHANDLER_H

// Includes
#include <stdlib.h>
#include <TokenizedLineReader.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * HeatOptionHandler handles the case where we want to use heat equation.
 */
class HeatOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	HeatOptionHandler() :
			OptionHandler("heat",
					"heat <value1> <value2>            "
							"The heat flux (in W nm-2) at the surface and the temperature in the bulk (Kelvin).\n") {
	}

	/**
	 * The destructor
	 */
	~HeatOptionHandler() {
	}

	/**
	 * This method will set the IOptions heatFlag, constTemperature, and bulkTemp
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The value for the temperature.
	 */
	bool handler(IOptions *opt, const std::string& arg) {

		// Set the flag to use constant temperature to true
		opt->setHeatFlag(true);

		// Build an input stream from the argument string.
		xolotlCore::TokenizedLineReader<std::string> reader;
		auto argSS = std::make_shared < std::istringstream > (arg);
		reader.setInputStream(argSS);

		// Break the argument into tokens.
		auto tokens = reader.loadLine();

		// Set the value for the heat flux at the surface
		double temp = strtod(tokens[0].c_str(), NULL);
		opt->setConstTemperature(temp);

		// Set the value for the temperature in the bulk
		temp = strtod(tokens[1].c_str(), NULL);
		opt->setBulkTemperature(temp);

		return true;
	}

};
//end class HeatOptionHandler

} /* namespace xolotlCore */

#endif
