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
    			"material <material>               "
    			"This option allows the user to change the profile of "
    			"the helium flux corresponding to the material.  \n"
    			"                                    The material options are as follows: "
    			"{W100, W110, W111}, where W is for \n"
    			"                                    tungsten and the numbers correspond to "
    			"the surface orientation.\n")
	{}

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
    	// Set the material name
    	opt->setMaterial(arg);
    	return true;
    }

};//end class MaterialOptionHandler

} /* namespace xolotlCore */

#endif
