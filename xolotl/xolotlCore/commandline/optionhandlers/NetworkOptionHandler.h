#ifndef NETWORKOPTIONHANDLER_H
#define NETWORKOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * NetworkOptionHandler handles the name of the network file.
 */
class NetworkOptionHandler : public OptionHandler {
public:

	/**
	 * The default constructor
	 */
    NetworkOptionHandler() :
    	OptionHandler("networkFile",
    			"networkFile <filename>      "
    			"The network will be loaded from this HDF5 file.") {}

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
