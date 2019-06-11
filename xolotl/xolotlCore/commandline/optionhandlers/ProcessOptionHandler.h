#ifndef PROCESSOPTIONHANDLER_H
#define PROCESSOPTIONHANDLER_H

// Includes
#include <TokenizedLineReader.h>
#include <memory>
#include <string.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * ProcessOptionHandler handles the choice of physical processes for the simulation.
 */
class ProcessOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	ProcessOptionHandler() :
			OptionHandler("process",
					"process                           "
							"List of all the processes to use in the simulation (reaction, diff, advec, "
							"modifiedTM, movingSurface, "
					"\n	                            bursting, attenuation, resolution).\n") {
	}

	/**
	 * The destructor
	 */
	~ProcessOptionHandler() {
	}

	/**
	 * This method will set the IOptions processMap to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The list of processes.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Initialize the map of processes
		std::map<std::string, bool> processMap;
		processMap["reaction"] = false;
		processMap["diff"] = false;
		processMap["advec"] = false;
		processMap["modifiedTM"] = false;
		processMap["movingSurface"] = false;
		processMap["bursting"] = false;
		processMap["attenuation"] = false;
		processMap["resolution"] = false;

		// Build an input stream from the argument string.
		xolotlCore::TokenizedLineReader<std::string> reader;
		auto argSS = std::make_shared<std::istringstream>(arg);
		reader.setInputStream(argSS);

		// Break the argument into tokens.
		auto tokens = reader.loadLine();

		// Loop on the tokens
		for (int i = 0; i < tokens.size(); ++i) {
			// Switch the value to true in the map
			processMap[tokens[i]] = true;
		}

		// Give it to the options
		opt->setProcesses(processMap);

		return true;
	}

};
//end class ProcessOptionHandler

} /* namespace xolotlCore */

#endif
