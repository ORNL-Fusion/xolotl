#ifndef PETSCOPTIONHANDLER_H
#define PETSCOPTIONHANDLER_H

// Includes
#include <TokenizedLineReader.h>
#include <memory>
#include <string.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * PetscOptionHandler handles all the Petsc options.
 */
class PetscOptionHandler: public OptionHandler {
protected:

	/**
	 * The default constructor
	 */
	PetscOptionHandler() :
			OptionHandler() {
	}
	;

public:

	/**
	 * The constructor to use.
	 * @param keyName The name for the key.
	 * @param msg The help message.
	 */
	PetscOptionHandler(std::string keyName, std::string msg) :
			OptionHandler(keyName, msg) {
	}

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
	bool handler(IOptions *opt, std::string arg) {
		// Build an input stream from the argument string.
		xolotlCore::TokenizedLineReader<std::string> reader;
		auto argSS = std::make_shared < std::istringstream > (arg);
		reader.setInputStream(argSS);

		// Break the argument into tokens.
		auto tokens = reader.loadLine();

		// Construct the argv from the stream of tokens.
		opt->setPetscArgc(tokens.size());
		// The PETSc argv is an array of pointers to C strings.
		auto argv = new char*[tokens.size() + 1];
		int idx = 0;
		for (auto iter = tokens.begin(); iter != tokens.end(); ++iter) {
			argv[idx] = new char[iter->length() + 1];
			strcpy(argv[idx], iter->c_str());
			++idx;
		}
		argv[idx] = 0; // null-terminate the array

		// Set the petscArgv
		opt->setPetscArgv(argv);

		return true;
	}

};
//end class PetscOptionHandler

} /* namespace xolotlCore */

#endif
