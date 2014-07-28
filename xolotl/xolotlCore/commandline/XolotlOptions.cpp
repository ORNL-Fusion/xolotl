#include <iostream>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <string.h>
#include <TokenizedLineReader.h>
#include "XolotlOptions.h"

namespace xolotlCore {

XolotlOptions::XolotlOptions( void )
  : useConstTempHandlers ( false ),	// by default, use constant temperature handlers
    useTempProfileHandlers ( false ),
    usePerfStdHandlers( true ),  // by default, use "std" handlers for performance
    useVizStdHandlers( false ), // and "dummy" handlers for visualization
    useMaxHeFluence( false ),   // by default do not use the maximum Helium fluence option
    petscArgc( 0 ),
    petscArgv( NULL ),
    constTemp(1000)	// by default, the constant temperature is 1000K
{
    // Add our notion of which options we support.
    optionsMap["networkFile"] = new OptInfo(
        "networkFile                 The network will be loaded from this file.",
        handleNetworkOptionCB );
    optionsMap["startTemp"] = new OptInfo(
        "startTemp <const_temp>      The temperature (in Kelvin) will be the constant floating point value specified. (default = 1000)"
        "\n			      (NOTE: Use only ONE temperature option)",
        handleConstTemperatureOptionCB );
    optionsMap["tempFile"] = new OptInfo(
        "tempFile <filename>         A temperature profile is given by the specified file, then linear interpolation is used to fit the data."
        "\n			      (NOTE: If a temperature file is given, a constant temperature should NOT be given)",
        handleTemperatureFileOptionCB );
    optionsMap["perfHandler"] = new OptInfo(
        "perfHandler {std,dummy}     Which set of performance handlers to use. (default = std)",
        handlePerfHandlersOptionCB );
    optionsMap["vizHandler"] = new OptInfo(
        "vizHandler {std,dummy}      Which set of handlers to use for the visualization. (default = dummy)",
        handleVizHandlersOptionCB );
    optionsMap["maxHeFluence"] = new OptInfo(
        "maxHeFluence <value>      The maximum value of the Helium fluence the user wishes to integrate to.",
        handleHeFluenceOptionCB );
    optionsMap["petscArgs"] = new OptInfo(
        "petscArgs                   All the arguments that will be given to PETSc",
        handlePetscOptionCB );
}


XolotlOptions::~XolotlOptions( void )
{
    // release the dynamically-allocated PETSc arguments
    for( unsigned int i = 0; i < petscArgc; ++i )
    {
        delete[] petscArgv[i];
    }
    delete[] petscArgv;
    petscArgv = NULL;
}


void
XolotlOptions::readParams( int argc, char* argv[] )
{
	// Check if we were given at least our positional arguments.
    if( argc < 1 )
    {
        std::cerr << "Insufficient input provided! Aborting!" << std::endl;
        showHelp( std::cerr );
        shouldRunFlag = false;
        exitCode = EXIT_FAILURE;

        return;
    }
    else
    {
        // Try to interpret first argument as the parameters file name.
        // If it starts with a dash, it's probably not.
        // But it is a common thing to try to get help by passing
        // only a help flag, so we check for that special case.
        std::string currArg = argv[0];
        if( currArg == "--help" )
        {
        	// Print the help message asked by the user
        	showHelp( std::cout );
            shouldRunFlag = false;
        }
        else if( currArg[0] == '-' )
        {
            // this looks like an option, not the required filename.
            std::cerr << "No parameter file provided.  Aborting!" << std::endl;
            showHelp( std::cerr );
            shouldRunFlag = false;
            exitCode = EXIT_FAILURE;

            return;
        }
        else
        {
            // Let the base Options class handle our options.
            Options::readParams(argc, argv);
        }
    }

    return;
}


void
XolotlOptions::showHelp( std::ostream& os ) const
{
    os << "usage: xolotl param_file_name \n\n"
        << "See the Xolotl documentation for PETSc options. \n"
        << "param_file_name should contain:"
        << std::endl;
    Options::showHelp( os );
}


bool
XolotlOptions::handleNetworkOption( std::string arg )
{
    bool ret = true;

    // The base class should check for situations where
    // we expect an argument but don't get one.
    assert( !arg.empty() );

    // Set the name of the network file
    networkFileName = arg;

    return ret;
}

bool
XolotlOptions::handleNetworkOptionCB( Options* opts, std::string arg )
{
    return static_cast<XolotlOptions*>( opts )->handleNetworkOption( arg );
}


bool
XolotlOptions::handleConstTemperatureOption( std::string arg )
{
	bool ret = true;

    // The base class should check for situations where
    // we expect an argument but don't get one.
    assert( !arg.empty() );

    useConstTempHandlers = true;
    constTemp = strtod(arg.c_str(), NULL);
    tempProfileFileName = "";

    return ret;
}

bool
XolotlOptions::handleConstTemperatureOptionCB( Options* opts, std::string arg )
{
    return static_cast<XolotlOptions*>( opts )->handleConstTemperatureOption( arg );
}


bool
XolotlOptions::handleTemperatureFileOption( std::string arg )
{
    bool ret = true;

    // The base class should check for situations where
    // we expect an argument but don't get one.
    assert( !arg.empty() );

    // Check to make sure the temperature file exists
    std::ifstream inFile(arg.c_str());
    if (!inFile)
    {
    	std::cerr << "\nCould not open file containing temperature profile data.  "
    			"Aborting!" << std::endl;
    	showHelp( std::cerr );
    	shouldRunFlag = false;
    	exitCode = EXIT_FAILURE;
    	ret = false;
    }
    else
    {
        useTempProfileHandlers = true;
        tempProfileFileName = arg;
        constTemp = 0;
    }

    return ret;
}

bool
XolotlOptions::handleTemperatureFileOptionCB( Options* opts, std::string arg )
{
    return static_cast<XolotlOptions*>( opts )->handleTemperatureFileOption( arg );
}


bool
XolotlOptions::handlePerfHandlersOption( std::string arg )
{
    bool ret = true;

    // The base class should check for situations where
    // we expect an argument but don't get one.
    assert( !arg.empty() );

    // Determine the type of handlers we are being asked to use
    if( arg == "std" )
    {
        usePerfStdHandlers = true;
    }
    else if( arg == "dummy" )
    {
        usePerfStdHandlers = false;
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
XolotlOptions::handlePerfHandlersOptionCB( Options* opts, std::string arg )
{
    return static_cast<XolotlOptions*>( opts )->handlePerfHandlersOption( arg );
}


bool
XolotlOptions::handleVizHandlersOption( std::string arg )
{
    bool ret = true;

    // The base class should check for situations where
    // we expect an argument but don't get one.
    assert( !arg.empty() );

    // Determine the type of handlers we are being asked to use
    if( arg == "std" )
    {
        useVizStdHandlers = true;
    }
    else if( arg == "dummy" )
    {
        useVizStdHandlers = false;
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
XolotlOptions::handleVizHandlersOptionCB( Options* opts, std::string arg )
{
    return static_cast<XolotlOptions*>( opts )->handleVizHandlersOption( arg );
}

bool
XolotlOptions::handleHeFluenceOption( std::string arg )
{
    bool ret = true;

    // The base class should check for situations where
    // we expect an argument but don't get one.
    assert( !arg.empty() );

    useMaxHeFluence = true;
    maxHeliumFluence = strtod(arg.c_str(), NULL);

    return ret;
}

bool
XolotlOptions::handleHeFluenceOptionCB( Options* opts, std::string arg )
{
    return static_cast<XolotlOptions*>( opts )->handleHeFluenceOption( arg );
}


bool
XolotlOptions::handlePetscOption( std::string arg )
{
    bool ret = true;
    assert( !arg.empty() );

    // Build an input stream from the argument string.
    xolotlCore::TokenizedLineReader<std::string> reader;
    auto argSS = std::make_shared<std::istringstream>(arg);
    reader.setInputStream(argSS);

    // Break the argument into tokens.
    auto tokens = reader.loadLine();

    // Construct the PETSc argv from the stream of tokens.
    // The PETSc argv is an array of pointers to C strings.
    petscArgc = tokens.size();
    petscArgv = new char*[petscArgc+1];
    int idx = 0;
    for( auto iter = tokens.begin(); iter != tokens.end(); ++iter )
    {
        petscArgv[idx] = new char[iter->length()+1];
        strcpy(petscArgv[idx], iter->c_str());
        ++idx;
    }
    petscArgv[idx] = 0; // null-terminate the array

    return ret;
}


bool
XolotlOptions::handlePetscOptionCB( Options* opts, std::string arg )
{
    return static_cast<XolotlOptions*>( opts )->handlePetscOption( arg );
}

}; // end namespace xolotlCore


