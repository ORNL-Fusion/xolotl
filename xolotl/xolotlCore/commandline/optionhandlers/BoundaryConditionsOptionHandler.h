#ifndef BOUNDARYCONDITIONSOPTIONHANDLER_H
#define BOUNDARYCONDITIONSOPTIONHANDLER_H

// Includes
#include <string.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * BoundaryConditionsOptionHandler handles the boundary conditions on the left and right side of the grid.
 */
class BoundaryConditionsOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	BoundaryConditionsOptionHandler() :
			OptionHandler("boundary",
					"boundary <values>                 "
							"This option allows the user to choose the boundary conditions in the X direction.\n"
							"                                    The first one correspond to the left side (surface) and second one to the right (bulk), "
							"0 means mirror, 1 means free surface.\n") {
	}

	/**
	 * The destructor
	 */
	~BoundaryConditionsOptionHandler() {
	}

	/**
	 * This method will set the IOptions the left and right boundaries
	 * to the value given as the argument.
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

		// Set the left boundary
		opt->setLeftBoundary(strtol(tokens[0].c_str(), NULL, 10));

		// Set the right boundary
		opt->setRightBoundary(strtol(tokens[1].c_str(), NULL, 10));

		return true;
	}

};
//end class BoundaryConditionsOptionHandler

} /* namespace xolotlCore */

#endif
