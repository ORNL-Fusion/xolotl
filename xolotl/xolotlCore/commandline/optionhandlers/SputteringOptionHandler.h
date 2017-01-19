#ifndef SPUTTERINGOPTIONHANDLER_H
#define SPUTTERINGOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * SputteringOptionHandler handles the number of sputtering yield option.
 */
class SputteringOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	SputteringOptionHandler() :
			OptionHandler("sputtering",
					"sputtering <yield>                "
							"This option allows the user to add a sputtering yield (atoms/ion)") {
	}

	/**
	 * The destructor
	 */
	~SputteringOptionHandler() {
	}

	/**
	 * This method will set the IOptions sputteringYield
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The sputtering yield.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Convert to double
		double yield = strtod(arg.c_str(), NULL);;
		// Set the sputtering yield
		opt->setSputteringYield(yield);

		return true;
	}

};
//end class SputteringOptionHandler

} /* namespace xolotlCore */

#endif
