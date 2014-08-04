#ifndef CONSTTEMPOPTIONHANDLER_H
#define CONSTTEMPOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * ConstTempOptionHandler handles the case where we want to use a constant temperature.
 */
class ConstTempOptionHandler : public OptionHandler {
protected:

	/**
	 * The default constructor
	 */
    ConstTempOptionHandler() : OptionHandler() {};


public:

	/**
	 * The constructor to use.
	 * @param keyName The name for the key.
	 * @param msg The help message.
	 */
    ConstTempOptionHandler(std::string keyName, std::string msg) :
    	OptionHandler(keyName, msg) {}

	/**
	 * The destructor
	 */
    ~ConstTempOptionHandler() {}

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

};//end class ConstTempOptionHandler

} /* namespace xolotlCore */

#endif
