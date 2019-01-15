#ifndef ESTOPPINGPOWEROPTIONHANDLER_H
#define ESTOPPINGPOWEROPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * EStoppingPowerOptionHandler handles the choice of the re-solution fit in
 * the material option.
 */
class EStoppingPowerOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	EStoppingPowerOptionHandler() :
			OptionHandler("zeta",
					"zeta <value>                      "
							"The value of the electronic stopping power in the material (0.73 by default).\n") {
	}

	/**
	 * The destructor
	 */
	~EStoppingPowerOptionHandler() {
	}

	/**
	 * This method will set the IOptions Zeta
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The value for the initial vacancy concentration.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Set the value for the initial vacancy concentration
		double zeta = strtod(arg.c_str(), NULL);

		opt->setZeta(zeta);
		return true;
	}

};
//end class EStoppingPowerOptionHandler

} /* namespace xolotlCore */

#endif
