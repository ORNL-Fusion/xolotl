#ifndef FLUXOPTIONHANDLER_H
#define FLUXOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * FluxOptionHandler handles the flux option.
 */
class FluxOptionHandler : public OptionHandler {
protected:

	/**
	 * The default constructor
	 */
    FluxOptionHandler() : OptionHandler() {};


public:

	/**
	 * The constructor to use.
	 * @param keyName The name for the key.
	 * @param msg The help message.
	 */
    FluxOptionHandler(std::string keyName, std::string msg) :
    	OptionHandler(keyName, msg) {}

	/**
	 * The destructor
	 */
    ~FluxOptionHandler() {}

    /**
     * This method will set the IOptions heliumFluenceFlag and maxHeliumFluence
     * to the value given as the argument.
     *
     * @param opt The pointer to the option that will be modified.
     * @param arg The value for the maximum fluence.
     */
    bool handler(IOptions *opt, std::string arg) {
    	// Set the corresponding flag to true
    	opt->setHeliumFluxFlag(true);

    	// Set the value for the maximum fluence
    	double flux = strtod(arg.c_str(), NULL);

    	opt->setHeliumFlux(flux);
    	return true;
    }

};//end class FluxOptionHandler

} /* namespace xolotlCore */

#endif
