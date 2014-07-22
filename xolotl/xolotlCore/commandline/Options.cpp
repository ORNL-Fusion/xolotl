#include <string>
#include <cassert>
#include <memory>
#include <iostream>
#include <fstream>
#include "Options.h"


namespace xolotlCore {


Options::Options( void )
  : shouldRunFlag( true ),
    exitCode( EXIT_SUCCESS )
{
}


Options::~Options( void )
{
    // Release the items in our map of potential options.
    for( auto iter = optionsMap.begin(); iter != optionsMap.end(); iter++ )
    {
        OptInfo* currOptInfo = iter->second;
        delete currOptInfo;
    }
    optionsMap.clear();
}



void
Options::showHelp( std::ostream& os ) const
{
    os << "Supported options:\n";

    for( OptionsMap::const_iterator iter = optionsMap.begin(); iter != optionsMap.end(); iter++ )
    {
        os << "  " << iter->second->helpMessage << '\n';
    }
    os << std::endl;
}

std::vector<std::string>
Options::ParamReadLine(std::shared_ptr<std::istream> inputstream) {
	// Local Declarations
	size_t lastDelimiterPos = 0, nextDelimiterPos = 0, finalDelimiterPos =
			std::string::npos;
	std::vector<std::string> paramVector;
	std::string subLine;

	// Check the inputstream before reading from it
	if (!inputstream->eof() && inputstream->good()) {
		// Get the line
		std::string line;
		std::getline(*inputstream, line);
		// Split it if it is not empty and does not start with the comment
		// character
		if (!line.empty()) {
			// If this line is a comment, skip it by calling this operation
			// again
			if (line.find("#") == 0)
				return ParamReadLine(inputstream);
			// Remove delimiters at the beginning of the string
			if (line.find("=") == 0)
				line = line.substr(1);
			// Remove delimiters at the end of the string
			if (line.find("=", line.size() - 2)
					== line.size() - 1)
				line = line.erase(line.size() - 1);
			// Find the first instance of the delimiter
			nextDelimiterPos = line.find("=");
			// Only split the line if it contains the delimiter
			if (nextDelimiterPos != finalDelimiterPos) {
				// Walk across each piece of data in the line stopping only
				// when the end of the string is reached.
				while (lastDelimiterPos != finalDelimiterPos) {
					// Create the substring that holds the single piece of
					// data between the delimiters
					subLine = line.substr(lastDelimiterPos,
							nextDelimiterPos - lastDelimiterPos);
					paramVector.push_back(subLine);
					// Switch the delimiter positions and find the next
					lastDelimiterPos =
							(nextDelimiterPos != finalDelimiterPos) ?
									nextDelimiterPos + 1 : nextDelimiterPos;
					nextDelimiterPos = line.find("=",
							lastDelimiterPos);
				}
			} else {
				// Otherwise just put an empty line in the array
				paramVector.push_back("");
			}
		}
	}

	return paramVector;
}


void
Options::readParams( int argc, char* argv[] )
{
    // All the options are read from an ASCII file that is parsed
	// with the TokenizedLineReader.
	// We assume that the name of this file is the first and only
	// argument.

	// Load the content of the file in a stream

	// Create the param stream
	std::shared_ptr<std::ifstream> paramStream;
	paramStream = std::make_shared<std::ifstream>(argv[0]);

	if (!paramStream->good()) {
        // The file is empty.
        std::cerr << "The parameter file is empty. Aborting!" << std::endl;
        showHelp( std::cerr );
        shouldRunFlag = false;
        exitCode = EXIT_FAILURE;

        return;
	}

	// Local Declarations
	// Load the first line
	auto line = ParamReadLine(paramStream);
	// And start looping on the lines
	while (line.size() > 0) {
		auto iter = optionsMap.find(line[0]);
		// If the option if found
		if(iter != optionsMap.end()) {
            // Call the option's handler
            OptInfo* currOptInfo = iter->second;
            assert( currOptInfo != NULL );
            // Continue to read if everything went well with the current option
            bool continueReading = (*currOptInfo->optHandler)( this, line[1] );

            if( !continueReading )
            {
            	// Something went wrong.
                std::cerr << "Option: Something went wrong setting the options." << std::endl;
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

		line = ParamReadLine(paramStream);
	}

    return;
}

}; // end namespace xolotlCore


