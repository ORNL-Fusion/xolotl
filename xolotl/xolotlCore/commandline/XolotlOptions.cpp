#include <cassert>
#include "xolotlCore/commandline/XolotlOptions.h"

namespace xolotlCore {


XolotlOptions::XolotlOptions( void )
  : useStdHandlers( true )              // by default, use "standard" handlers
{
    // Add our notion of which options we support.
    optionsMap["--handlers"] = new OptInfo( true,
        "--handlers {std,dummy}     Which set of handlers to use.",
        handleHandlersOptionCB );
    optionsMap["--petsc"] = new OptInfo( false,
        "--petsc                    All subsequent command line args should be given to PETSc",
        handlePetscOptionCB );
}


int
XolotlOptions::parseCommandLine( int argc, char* argv[] )
{
    int nArgsUsed = 0;

    // Check if we were given at least our positional arguments.
    if( argc < 1 )
    {
        std::cerr << "Insufficient input provided! Aborting!" << std::endl;
        showHelp( std::cerr );
        shouldRunFlag = false;
        exitCode = EXIT_FAILURE;
    }
    else
    {
        // Try to interpret first argument as the network file name.
        // If it starts with a dash, it's probably not.
        // But it is a common thing to try to get help by passing
        // only a help flag, so we check for that special case.
        std::string currArg = argv[0];
        if( currArg == "--help" )
        {
            // let the base class handle it
            nArgsUsed = Options::parseCommandLine( argc, argv );
        }
        else if( currArg[0] == '-' )
        {
            // this looks like an option, not the required filename.
            std::cerr << "No network file provided.  Aborting!" << std::endl;
            showHelp( std::cerr );
            shouldRunFlag = false;
            exitCode = EXIT_FAILURE;
        }
        else
        {
            // Use the first argument as the network file name
            netFileName = argv[0];
            nArgsUsed = 1;  // one for network file name

            // Let the base Options class handle our options.
            nArgsUsed += Options::parseCommandLine( argc - nArgsUsed, argv + nArgsUsed );
        }
    }

    return nArgsUsed;
}


void
XolotlOptions::showHelp( std::ostream& os ) const
{
    os << "usage: xolotl network_file_name [OPTIONS]\n\n"
        << "See the Xolotl documentation for PETSc options."
        << std::endl;
    Options::showHelp( os );
}


bool
XolotlOptions::handleHandlersOption( std::string arg )
{
    bool ret = true;

    // The base class should check for situations where
    // we expect an argument but don't get one.
    assert( !arg.empty() );

    // Determine the type of handlers we are being asked to use
    if( arg == "std" )
    {
        useStdHandlers = true;
    }
    else if( arg == "dummy" )
    {
        useStdHandlers = false;
    }
    else
    {
        std::cerr << "Options: unrecognized argument " << arg << std::endl;
        showHelp( std::cerr );
        shouldRunFlag = false;
        exitCode = EXIT_FAILURE;
        ret = false;
    }

    return ret;
}

bool
XolotlOptions::handleHandlersOptionCB( Options* opts, std::string arg )
{
    return static_cast<XolotlOptions*>( opts )->handleHandlersOption( arg );
}


bool
XolotlOptions::handlePetscOption( std::string arg )
{
    assert( arg.empty() );

    // we are done parsing our own arguments, the rest are assumed 
    // to be arguments for PETSc.
    return false;   
}


bool
XolotlOptions::handlePetscOptionCB( Options* opts, std::string arg )
{
    return static_cast<XolotlOptions*>( opts )->handlePetscOption( arg );
}

}; // end namespace xolotlCore


