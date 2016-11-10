#ifndef OPTIONHANDLER_H
#define OPTIONHANDLER_H

// Includes
#include "IOptionHandler.h"

namespace xolotlCore {

/**
 * OptionHandler realizes the IOptionHandler interface.
 * There will be an implementation for each different option.
 */
class OptionHandler: public IOptionHandler {
protected:

	/**
	 * The default constructor
	 */
	OptionHandler();

public:

	/**
	 * The constructor to use.
	 * @param keyName The name for the key.
	 * @param msg The help message.
	 */
	OptionHandler(const std::string& keyName, const std::string& msg) {
		key = keyName;
		helpMessage = msg;
	}

	/**
	 * The destructor
	 */
	~OptionHandler() {
	}

	/**
	 * The function that will handle the specific option.
	 * Every subclass will have to implement this function.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The argument for the option.
	 */
	virtual bool handler(IOptions * opt, const std::string& arg) {
		return false;
	}

};
//end class OptionHandler

} /* namespace xolotlCore */

#endif
