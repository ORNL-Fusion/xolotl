#ifndef GRIDPARAMOPTIONHANDLER_H
#define GRIDPARAMOPTIONHANDLER_H

// Includes
#include <string.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * GridParamOptionHandler handles the option to generate
 * the grid if not loaded from the HDF5 file.
 */
class GridParamOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	GridParamOptionHandler() :
			OptionHandler("grid",
					"grid                              "
							"This option allows the user to define the boundaries of the grid.\n"
							"                                    To do so, simply write the values in order "
							"nX xStepSize nY yStepSize nZ zStepSize .\n") {
	}

	/**
	 * The destructor
	 */
	~GridParamOptionHandler() {
	}

	/**
	 * This method will set the IOptions grid parameters
	 * to the values given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The list of PETSc options.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Build an input stream from the argument string.
		xolotlCore::TokenizedLineReader<std::string> reader;
		auto argSS = std::make_shared<std::istringstream>(arg);
		reader.setInputStream(argSS);

		// Break the argument into tokens.
		auto tokens = reader.loadLine();

		// Set the values for the for the depth
		int n = strtol(tokens[0].c_str(), NULL, 10);
		opt->setNX(n);
		double size = strtod(tokens[1].c_str(), NULL);
		opt->setXStepSize(size);

		// Check if we have other values
		if (tokens.size() > 2) {
			// Set the values for the for Y
			n = strtol(tokens[2].c_str(), NULL, 10);
			opt->setNY(n);
			size = strtod(tokens[3].c_str(), NULL);
			opt->setYStepSize(size);

			// Check if we have other values
			if (tokens.size() > 4) {
				// Set the values for the for Z
				n = strtol(tokens[4].c_str(), NULL, 10);
				opt->setNZ(n);
				size = strtod(tokens[5].c_str(), NULL);
				opt->setZStepSize(size);
			}
		}

		return true;
	}

};
//end class GridParamOptionHandler

} /* namespace xolotlCore */

#endif
