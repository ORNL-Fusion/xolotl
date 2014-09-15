#ifndef PERFOPTIONHANDLER_H
#define PERFOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * PerfOptionHandler handles the choice of handlers for the performance infrastructure.
 */
class PerfOptionHandler: public OptionHandler {
public:

	/**
	 * Construct a PerfOptionHandler.
	 */
	PerfOptionHandler() :
		OptionHandler("perfHandler",
				"perfHandler {std,dummy,os,papi}     "
				"Which set of performance handlers to use. (default = std)") {}

	/**
	 * Destroy the PerfOptionHandler.
	 */
	~PerfOptionHandler() {
	}

	/**
	 * This method will set the IOptions perfStandardHandlersFlag
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The argument for the flag.
	 */
	bool handler(IOptions *opt, std::string arg) {
        
        bool ret = true;

        try
        {
		    // Determine the type of handlers we are being asked to use
            xolotlPerf::IHandlerRegistry::RegistryType rtype = xolotlPerf::toPerfRegistryType(arg);
            opt->setPerfHandlerType( rtype );
        }
        catch( std::invalid_argument& e )
        {
            std::cerr << e.what() << std::endl;
            opt->showHelp(std::cerr);
            opt->setShouldRunFlag(false);
            opt->setExitCode(EXIT_FAILURE);
            ret = false;
        }

		return ret;
	}

};
//end class PerfOptionHandler

} /* namespace xolotlCore */

#endif
