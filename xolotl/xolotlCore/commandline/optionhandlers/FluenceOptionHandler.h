#ifndef FLUENCEOPTIONHANDLER_H
#define FLUENCEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * FluenceOptionHandler handles the use of the maximum fluence option.
 */
class FluenceOptionHandler : public OptionHandler {
protected:

	/**
	 * The default constructor
	 */
    FluenceOptionHandler() : OptionHandler() {};


public:

	/**
	 * The constructor to use.
	 * @param keyName The name for the key.
	 * @param msg The help message.
	 */
    FluenceOptionHandler(std::string keyName, std::string msg) :
    	OptionHandler(keyName, msg) {}

	/**
	 * The destructor
	 */
    ~FluenceOptionHandler() {}

    /**
     * This method will set the IOptions heliumFluenceFlag and maxHeliumFluence
     * to the value given as the argument.
     *
     * @param opt The pointer to the option that will be modified.
     * @param arg The value for the maximum fluence.
     */
    bool handler(IOptions *opt, std::string arg) {
    	// Set the corresponding flag to true
    	opt->setHeliumFluenceFlag(true);

    	// Set the value for the maximum fluence
    	double fluence = strtod(arg.c_str(), NULL);

    	opt->setMaxHeliumFluence(fluence);
    	return true;
    }

};//end class FluenceOptionHandler

} /* namespace xolotlCore */

#endif
