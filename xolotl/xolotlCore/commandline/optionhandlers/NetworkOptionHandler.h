#ifndef NETWORKOPTIONHANDLER_H
#define NETWORKOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * NetworkOptionHandler handles the name of the network file.
 */
class NetworkOptionHandler : public OptionHandler {
protected:

	/**
	 * The default constructor
	 */
    NetworkOptionHandler() : OptionHandler() {};


public:

	/**
	 * The constructor to use.
	 * @param keyName The name for the key.
	 * @param msg The help message.
	 */
    NetworkOptionHandler(std::string keyName, std::string msg) :
    	OptionHandler(keyName, msg) {}

	/**
	 * The destructor
	 */
    ~NetworkOptionHandler() {}

    /**
     * This method will set the IOptions networkFilename
     * to the value given as the argument.
     *
     * @param opt The pointer to the option that will be modified.
     * @param arg The argument for the networkFilename.
     */
    bool handler(IOptions *opt, std::string arg) {
    	// Set the name of the network file
    	opt->setNetworkFilename(arg);

    	return true;
    }

};//end class NetworkOptionHandler

} /* namespace xolotlCore */

#endif
