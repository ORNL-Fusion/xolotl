#ifndef FLUXPROFILEOPTIONHANDLER_H
#define FLUXPROFILEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"
#include <mpi.h>

namespace xolotlCore {

/**
 * FluxProfileOptionHandler handles the case where we want to use a helium flux profile.
 */
class FluxProfileOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	FluxProfileOptionHandler() :
			OptionHandler("fluxFile",
					"fluxFile <filename>               "
							"A time profile for the flux is given by the specified file, "
							"then linear interpolation is used to fit the data."
							"\n	                            (NOTE: If a flux profile file is given, "
							"a constant helium flux should NOT be given)\n") {
	}

	/**
	 * The destructor
	 */
	~FluxProfileOptionHandler() {
	}

	/**
	 * This method will set the IOptions FluxProfileFlag and FluxProfileFilename
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The name of the file where the profile is stored.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		bool ret = true;

		// Check that the profile file exists
		std::ifstream inFile(arg.c_str());
		if (!inFile) {
			std::cerr
					<< "\nCould not open file containing time profile data for the flux. "
							"Aborting!\n" << std::endl;
			opt->setShouldRunFlag(false);
			opt->setExitCode(EXIT_FAILURE);
			ret = false;
		} else {
			// Set the flag to use a helium profile to true
			opt->setFluxProfileFlag(true);

			// Set the name of the file
			opt->setFluxProfileName(arg);
		}

		return ret;
	}

};
//end class FluxProfileOptionHandler

} /* namespace xolotlCore */

#endif
