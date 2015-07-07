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
		// Build an input stream from the argument string.
		xolotlCore::TokenizedLineReader<std::string> reader;
		auto argSS = std::make_shared < std::istringstream > (arg);
		reader.setInputStream(argSS);

		// Break the argument into tokens.
		auto tokens = reader.loadLine();

		// PETSc assumes that argv[0] in the arguments it is given is the
	    // program name.  But our parsing of the PETSc arguments from
	    // the input parameter file gives us only the PETSc arguments without
	    // the program name as argv[0].  So - we adjust the arguments array.

		// Construct the argv from the stream of tokens.
		opt->setPetscArgc(tokens.size() + 1);

		// The PETSc argv is an array of pointers to C strings.
		auto argv = new char*[tokens.size() + 2];
		// Create the fake application name
		std::string appName = "fakeXolotlApplicationNameForPETSc";
		argv[0] = new char[appName.length() + 1];
		strcpy(argv[0], appName.c_str());

		// Now loop on the actual PETSc options
		int idx = 1;
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
