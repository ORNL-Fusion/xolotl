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
							"This option allows the user to choose the boundary conditions.\n"
							"                                    The first one correspond to the left side (surface) "
							"and second one to the right (bulk),\n"
							"                                    then two for Y and two for Z. "
							"0 means mirror or periodic, 1 means free surface.\n") {
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

		if (tokens.size() > 2)
			// Set the bottom boundary
			opt->setBottomBoundary(strtol(tokens[2].c_str(), NULL, 10));

		if (tokens.size() > 3)
			// Set the top boundary
			opt->setTopBoundary(strtol(tokens[3].c_str(), NULL, 10));

		if (tokens.size() > 4)
			// Set the front boundary
			opt->setFrontBoundary(strtol(tokens[4].c_str(), NULL, 10));

		if (tokens.size() > 5)
			// Set the back boundary
			opt->setBackBoundary(strtol(tokens[5].c_str(), NULL, 10));

		return true;
	}

};
//end class BoundaryConditionsOptionHandler

} /* namespace xolotlCore */

#endif
