#ifndef NETWORKPARAMOPTIONHANDLER_H
#define NETWORKPARAMOPTIONHANDLER_H

// Includes
#include <string.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * NetworkParamOptionHandler handles the option to generate
 * the network if not loaded from the HDF5 file.
 */
class NetworkParamOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	NetworkParamOptionHandler() :
			OptionHandler("netParam",
					"netParam                          "
							"This option allows the user to define the boundaries of the network.\n"
							"                                    To do so, simply write the values in order "
							"maxHe/Xe maxD maxT maxV maxI bool .\n") {
	}

	/**
	 * The destructor
	 */
	~NetworkParamOptionHandler() {
	}

	/**
	 * This method will set the IOptions network parameters
	 * to the values given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The list of PETSc options.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Set the flag to not use the HDF5 file
		opt->setHDF5Flag(false);

		// Build an input stream from the argument string.
		xolotlCore::TokenizedLineReader<std::string> reader;
		auto argSS = std::make_shared<std::istringstream>(arg);
		reader.setInputStream(argSS);

		// Break the argument into tokens.
		auto tokens = reader.loadLine();

		// Set the value for the impurities
		int max = strtol(tokens[0].c_str(), NULL, 10);
		opt->setMaxImpurity(max);

		// Check if we have other values
		if (tokens.size() > 1) {
			// Set the deuterium size
			max = strtol(tokens[1].c_str(), NULL, 10);
			opt->setMaxD(max);
			// Set the tritium size
			max = strtol(tokens[2].c_str(), NULL, 10);
			opt->setMaxT(max);
			// Set the vacancy size
			max = strtol(tokens[3].c_str(), NULL, 10);
			opt->setMaxV(max);
			// Set the interstitial size
			max = strtol(tokens[4].c_str(), NULL, 10);
			opt->setMaxI(max);

			// Check if there are other values
			if (tokens.size() > 5) {
				// Set the phase cut
				opt->setPhaseCutFlag(tokens[5] == "true");
			}
		}

		return true;
	}

}
;
//end class NetworkParamOptionHandler

} /* namespace xolotlCore */

#endif
