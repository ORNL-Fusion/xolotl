#ifndef VCONCENTRATIONOPTIONHANDLER_H
#define VCONCENTRATIONOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * VConcentrationOptionHandler handles the use of the initial vacancy concentration in
 * the material option.
 */
class VConcentrationOptionHandler : public OptionHandler {
public:

	/**
	 * The default constructor
	 */
    VConcentrationOptionHandler() :
    	OptionHandler("initialV",
    			"initialV <value>        "
    			"The value of the initial concentration of vacancies in the material.") {}

	/**
	 * The destructor
	 */
    ~VConcentrationOptionHandler() {}

    /**
     * This method will set the IOptions initialVConcentration
     * to the value given as the argument.
     *
     * @param opt The pointer to the option that will be modified.
     * @param arg The value for the initial vacancy concentration.
     */
    bool handler(IOptions *opt, std::string arg) {
    	// Set the value for the initial vacancy concentration
    	double conc = strtod(arg.c_str(), NULL);

    	opt->setInitialVConcentration(conc);
    	return true;
    }

};//end class VConcentrationOptionHandler

} /* namespace xolotlCore */

#endif
