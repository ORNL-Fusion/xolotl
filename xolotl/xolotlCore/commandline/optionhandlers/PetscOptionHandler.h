#ifndef PETSCOPTIONHANDLER_H
#define PETSCOPTIONHANDLER_H

// Includes
#include <TokenizedLineReader.h>
#include <memory>
#include <string.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * PetscOptionHandler handles all the PETSc options.
 */
class PetscOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	PetscOptionHandler() :
		OptionHandler("petscArgs",
				"petscArgs                         "
				"All the arguments that will be given to PETSc.\n") {}

	/**
	 * The destructor
	 */
	~PetscOptionHandler() {
	}

	/**
	 * This method will set the IOptions petscArgc and petscArgv
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The list of PETSc options.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Set the petscArgv that is the string of options
		opt->setPetscArgv(arg);

		return true;
	}

};
//end class PetscOptionHandler

} /* namespace xolotlCore */

#endif
