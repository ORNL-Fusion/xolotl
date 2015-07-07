#ifndef IOPTIONHANDLER_H
#define IOPTIONHANDLER_H

// Includes
#include <string>
#include <IOptions.h>

namespace xolotlCore {

/**
 * IOptionHandler describes the structure needed for the option handlers.
 * There will be an implementation for each different option.
 */
class IOptionHandler {
protected:

public:

	/**
	 * The description of how to use this option.
	 */
	std::string helpMessage;

	/**
	 * The name of the option as it would appear in the param.txt file
	 */
	std::string key;

	/**
	 * The destructor
	 */
    virtual ~IOptionHandler() {}

    /**
     * The function that will handle the specific option.
     * Every subclass will have to implement this function.
     *
     * @param opt The pointer to the option that will be modified.
     * @param arg The argument for the option.
     */
    virtual bool handler(IOptions *opt, const std::string& arg) = 0;

};//end class IOptionHandler

} /* namespace xolotlCore */

#endif
