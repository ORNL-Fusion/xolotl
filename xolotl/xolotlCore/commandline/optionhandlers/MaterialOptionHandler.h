#ifndef MATERIALOPTIONHANDLER_H
#define MATERIALOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * MaterialOptionHandler handles the material option.
 */
class MaterialOptionHandler : public OptionHandler {
public:

	/**
	 * The default constructor
	 */
    MaterialOptionHandler() :
    	OptionHandler("material",
    			"material <material>         "
    			"This option allows the user to change the profile of "
    			"the helium flux corresponding to the material.  (default = W100) \n"
    			"                              The material options are as follows: "
    			"{W100, W110, W111, W210, W211, W221, W310, W311, W320, W321, Fe}") {}

	/**
	 * The destructor
	 */
    ~MaterialOptionHandler() {}

    /**
     * This method will set the IOptions materialFlag and materialName
     * to the value given as the argument.
     *
     * @param opt The pointer to the option that will be modified.
     * @param arg The name of the material.
     */
    bool handler(IOptions *opt, std::string arg) {
    	// Set the corresponding flag to true
    	opt->setMaterialFlag(true);

    	opt->setMaterial(arg);
    	return true;
    }

};//end class MaterialOptionHandler

} /* namespace xolotlCore */

#endif
