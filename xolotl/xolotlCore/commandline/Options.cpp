#include <string>
#include <cassert>
#include "xolotlCore/commandline/Options.h"


namespace xolotlCore {


Options::Options( void )
  : shouldRunFlag( true ),
    exitCode( EXIT_SUCCESS )
{
    // Construct our notion of which options we support.
    optionsMap["--help"] = new OptInfo( false,
        "--help                     Show this help message.",
        Options::handleHelpOptionCB );
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


int
Options::parseCommandLine( int argc, char* argv[] )
{
    // We are not using any command line parsing library,
    // so our parsing approach is a fairly brute force one.
    // 
    // We expect the command line to contain lots of options
    // that do not pertain to us (e.g., PETSc options).  But
    // we expect that those options will come after an 
    // option that tells us our options are done.
    unsigned int i = 0;
    while( i < argc )
    {
        std::string currArg = argv[i++];

        auto iter = optionsMap.find( currArg );
        if( iter != optionsMap.end() )
        {
            // We recognized the option.
            // Call the option's handler.
            OptInfo* currOptInfo = iter->second;
            assert( currOptInfo != NULL );

            std::string optArg;
            if( currOptInfo->argRequired )
            {
                if( i < argc )
                {
                    optArg = argv[i++];
                }
                else
                {
                    std::cerr << "Option: required argument missing for " << currArg << std::endl;
                    shouldRunFlag = false;
                    exitCode = EXIT_FAILURE;
                    break;
                }
            }
            bool continueParsing = (*currOptInfo->optHandler)( this, optArg );
            if( !continueParsing )
            {
                break;
            }
        }
        else
        {
            // We did not recognize the option.
            std::cerr << "Option: unrecognized option " << currArg << std::endl;
            shouldRunFlag = false;
            exitCode = EXIT_FAILURE;
            break;
        }
    }

    return i;
}


//
// Methods for handling specific options we support.
//

bool
Options::handleHelpOption( std::string arg )
{
    // Show the help message to standard output (since it
    // was requested via command line switch).
    showHelp( std::cout );

    // Indicate that the program shouldn't try to run.
    shouldRunFlag = false;

    // Stop parsing.
    return false;
}


bool
Options::handleHelpOptionCB( Options* opts, std::string arg )
{
    return opts->handleHelpOption( arg );
}

}; // end namespace xolotlCore


