#ifndef RNGOPTIONHANDLER_H
#define RNGOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * RNGOptionHandler handles the configuration of the random number generator
 * Xolotl will use.
 */
class RNGOptionHandler: public OptionHandler {
public:

	/**
	 * Construct a RNGOptionHandler.
	 */
	RNGOptionHandler() :
		OptionHandler("rng",
				"rng [print] [seed_value]   "
				"Allows user to specify seed used to initialize random number\n"
                "generator (default = determined from current time) and\n"
                "whether each process should print the seed value\n"
                "it uses (default = don't print)\n")
    {}

	/**
	 * Destroy the RNGOptionHandler.
	 */
	~RNGOptionHandler() {
	}

	/**
	 * This method will set the IOptions perfStandardHandlersFlag
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The argument for the flag.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
        
        bool ret = true;

        try
        {
            // Break string into tokens to be handled individually.
            xolotlCore::TokenizedLineReader<std::string> reader;
		    auto istr = std::make_shared<std::istringstream>(arg);
            reader.setInputStream(istr);
            auto tokens = reader.loadLine();
            size_t currIdx = 0;

            // Determine whether we should print the seed value.
            bool shouldPrintSeed = false;
            if(tokens[currIdx] == "print") {
                shouldPrintSeed = true;
                ++currIdx;
            }
            opt->setPrintRNGSeed(shouldPrintSeed);

            if(currIdx < tokens.size()) {
                // Convert arg to an integer.
                char* ep = NULL;
                auto useed = strtoul(tokens[currIdx].c_str(), &ep, 10);
                if(ep != (tokens[currIdx].c_str() + tokens[currIdx].length()))
                {
                    std::ostringstream estr;
                    estr << "Invalid random number generator seed \"" << arg << "\".  Must be a non-negative integer.";
                    throw std::invalid_argument(estr.str());
                }
                opt->setRNGSeed(useed);
            }
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << e.what() << std::endl;
            opt->showHelp(std::cerr);
            opt->setShouldRunFlag(false);
            opt->setExitCode(EXIT_FAILURE);
            ret = false;
        }

		return ret;
	}

};
//end class RNGOptionHandler

} /* namespace xolotlCore */

#endif
