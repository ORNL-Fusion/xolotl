#include <cassert>
#include <fstream>
#include <TokenizedLineReader.h>
#include <NetworkOptionHandler.h>
#include <PetscOptionHandler.h>
#include <ConstTempOptionHandler.h>
#include <TempProfileOptionHandler.h>
#include <FluenceOptionHandler.h>
#include <FluxOptionHandler.h>
#include <PerfOptionHandler.h>
#include <VizOptionHandler.h>
#include "Options.h"

namespace xolotlCore {

Options::Options() :
		shouldRunFlag(true),
		exitCode(EXIT_SUCCESS),
		petscArgc(0),
		petscArgv(NULL),
		constTempFlag(true),
		tempProfileFlag(false),
		constTemperature(1000.0),
		heliumFluenceFlag(false),
		heliumFluxFlag(false),
		perfStandardHandlersFlag(true),
		vizStandardHandlersFlag(false) {
	// Create the network option handler
	auto networkHandler = new
			NetworkOptionHandler("networkFile",
					"networkFile <filename>      The network will be loaded from this file.");
	// Create the PETSc option handler
	auto petscHandler = new
			PetscOptionHandler("petscArgs",
					"petscArgs                   All the arguments that will be given to PETSc");
	// Create the constant temperature option handler
	auto constTempHandler = new
			ConstTempOptionHandler("startTemp",
					"startTemp <value>           The temperature (in Kelvin) will be the constant floating point value specified. (default = 1000)"
							"\n	                      (NOTE: Use only ONE temperature option)");
	// Create the temperature profile option handler
	auto tempProfileHandler = new
			TempProfileOptionHandler("tempFile",
					"tempFile <filename>         A temperature profile is given by the specified file, then linear interpolation is used to fit the data."
							"\n	                      (NOTE: If a temperature file is given, a constant temperature should NOT be given)");
	// Create the maximum fluence option handler
	auto fluenceHandler = new
			FluenceOptionHandler("maxHeFluence",
					"maxHeFluence <value>        The maximum value of the Helium fluence the user wishes to integrate to.");
	// Create the flux option handler
	auto fluxHandler = new
			FluxOptionHandler("heFlux",
					"heFlux <value>              This option allows the user to change the Helium flux by the factor specified (in nm).");
	// Create the performance handler option handler
	auto perfHandler = new
			PerfOptionHandler("perfHandler",
					"perfHandler {std,dummy}     Which set of performance handlers to use. (default = std)");
	// Create the visualization handler option handler
	auto vizHandler = new
			VizOptionHandler("vizHandler",
					"vizHandler {std,dummy}      Which set of handlers to use for the visualization. (default = dummy)");

	// Add our notion of which options we support.
	optionsMap[networkHandler->key] = networkHandler;
	optionsMap[petscHandler->key] = petscHandler;
	optionsMap[constTempHandler->key] = constTempHandler;
	optionsMap[tempProfileHandler->key] = tempProfileHandler;
	optionsMap[fluenceHandler->key] = fluenceHandler;
	optionsMap[fluxHandler->key] = fluxHandler;
	optionsMap[perfHandler->key] = perfHandler;
	optionsMap[vizHandler->key] = vizHandler;
}

Options::~Options(void) {
	// Release the items in our map of potential options.
	for (auto iter = optionsMap.begin(); iter != optionsMap.end(); iter++) {
		auto currOpt = iter->second;
		delete currOpt;
	}
	optionsMap.clear();

	// release the dynamically-allocated PETSc arguments
	for (unsigned int i = 0; i < petscArgc; ++i) {
		delete[] petscArgv[i];
	}
	delete[] petscArgv;
	petscArgv = NULL;
}

void Options::readParams(int argc, char* argv[]) {
	// All the options are read from an ASCII file that is parsed
	// with the TokenizedLineReader.
	// We assume that the name of this file is the first and only
	// argument.

	// Load the content of the file in a stream
	// Create the param stream
	std::shared_ptr<std::ifstream> paramStream;
	paramStream = std::make_shared < std::ifstream > (argv[0]);

	if (!paramStream->good()) {
		// The file is empty.
		std::cerr << "The parameter file is empty. Aborting!" << std::endl;
		showHelp(std::cerr);
		shouldRunFlag = false;
		exitCode = EXIT_FAILURE;

		return;
	}

	// Declare the TokenizedLineReader to read the lines
	xolotlCore::TokenizedLineReader<std::string> reader;
	reader.setInputStream(paramStream);
	// Set the delimiter to "="
	reader.setDelimiter("=");

	// Load the first line
	auto line = reader.loadLine();
	// And start looping on the lines
	while (line.size() > 0) {
		auto iter = optionsMap.find(line[0]);
		// If the option if found
		if (iter != optionsMap.end()) {
			// Call the option's handler
			auto currOpt = iter->second;
			assert(currOpt != NULL);
			// Continue to read if everything went well with the current option
			bool continueReading = currOpt->handler(this, line[1]);

			if (!continueReading) {
				// Something went wrong.
				std::cerr << "Option: Something went wrong setting the options."
						<< std::endl;
				shouldRunFlag = false;
				exitCode = EXIT_FAILURE;
				break;
			}
		}

		else {
			// We did not recognize the option.
			std::cerr << "Option: unrecognized option " << line[0] << std::endl;
			shouldRunFlag = false;
			exitCode = EXIT_FAILURE;
			break;
		}

		line = reader.loadLine();
	}

	return;
}

void Options::showHelp(std::ostream& os) const {
	os << "usage: xolotl param_file_name \n\n"
			<< "See the Xolotl documentation for PETSc options. \n"
			<< "Supported options:\n";

	for (OptionsMap::const_iterator iter = optionsMap.begin();
			iter != optionsMap.end(); iter++) {
		os << "  " << iter->second->helpMessage << '\n';
	}
	os << std::endl;
}

}
;
// end namespace xolotlCore

