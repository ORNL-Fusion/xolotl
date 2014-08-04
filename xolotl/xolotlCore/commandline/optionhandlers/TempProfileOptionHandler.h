#ifndef TEMPPROFILEOPTIONHANDLER_H
#define TEMPPROFILEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * TempProfileOptionHandler handles the case where we want to use a temperature profile.
 */
class TempProfileOptionHandler: public OptionHandler {
protected:

	/**
	 * The default constructor
	 */
	TempProfileOptionHandler() :
			OptionHandler() {
	}
	;

public:

	/**
	 * The constructor to use.
	 * @param keyName The name for the key.
	 * @param msg The help message.
	 */
	TempProfileOptionHandler(std::string keyName, std::string msg) :
			OptionHandler(keyName, msg) {
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
		// Check the profile file exists
		std::ifstream inFile(arg.c_str());
		if (!inFile) {
			std::cerr
					<< "\nCould not open file containing temperature profile data. "
							"Aborting!" << std::endl;
			opt->showHelp(std::cerr);
			opt->setShouldRunFlag(false);
			opt->setExitCode(EXIT_FAILURE);

			return false;
		}

		// Set the flag to use a temperature profile to true
		opt->setTempProfileFlag(true);

		// Set the name of the file
		opt->setTempProfileFilename(arg);

		return true;
	}

};
//end class TempProfileOptionHandler

} /* namespace xolotlCore */

#endif
