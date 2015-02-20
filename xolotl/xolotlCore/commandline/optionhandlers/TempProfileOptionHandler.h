#ifndef TEMPPROFILEOPTIONHANDLER_H
#define TEMPPROFILEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"
#include <mpi.h>

namespace xolotlCore {

/**
 * TempProfileOptionHandler handles the case where we want to use a temperature profile.
 */
class TempProfileOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	TempProfileOptionHandler() :
			OptionHandler("tempFile",
					"tempFile <filename>               "
							"A temperature profile is given by the specified file, "
							"then linear interpolation is used to fit the data."
							"\n	                            (NOTE: If a temperature file is given, "
							"a constant temperature should NOT be given)\n") {
	}

	/**
	 * The destructor
	 */
	~TempProfileOptionHandler() {
	}

	/**
	 * This method will set the IOptions tempProfileFlag and tempProfileFilename
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The name of the file where the profile is stored.
	 */
	bool handler(IOptions *opt, std::string arg) {
		bool ret = true;

		// Check that the profile file exists
		std::ifstream inFile(arg.c_str());
		if (!inFile) {
			std::cerr
					<< "\nCould not open file containing temperature profile data. "
							"Aborting!\n" << std::endl;
			opt->setShouldRunFlag(false);
			opt->setExitCode(EXIT_FAILURE);
			ret = false;
		} else {
			// Set the flag to use a temperature profile to true
			opt->setTempProfileFlag(true);

			// Set the name of the file
			opt->setTempProfileFilename(arg);
		}

		return ret;
	}

};
//end class TempProfileOptionHandler

} /* namespace xolotlCore */

#endif
