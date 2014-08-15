#ifndef FLUENCEOPTIONHANDLER_H
#define FLUENCEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * FluenceOptionHandler handles the use of the maximum fluence option.
 */
class FluenceOptionHandler : public OptionHandler {
public:

	/**
	 * The default constructor
	 */
    FluenceOptionHandler() :
    	OptionHandler("maxHeFluence",
    			"maxHeFluence <value>        The maximum value of the Helium fluence the user wishes to integrate to.") {}

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
