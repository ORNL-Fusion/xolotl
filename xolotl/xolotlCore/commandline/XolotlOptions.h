#ifndef XOLOTLOPTIONS_H
#define XOLOTLOPTIONS_H

#include <string>
#include "Options.h"
#include "xolotlPerf/xolotlPerf.h"

namespace xolotlCore {

class XolotlOptions : public Options
{
private:

    // Use the constant temperature set of handlers?
    bool useConstTempHandlers;

    // Use the temperature profile set of handlers?
    bool useTempProfileHandlers;

    // Which type of performance infrastucture should we use?
    xolotlPerf::IHandlerRegistry::RegistryType perfRegistryType;

    // Use the "standard" set of handlers for the visualization?
    bool useVizStdHandlers;

    // Name of the input network file.
    std::string networkFileName;

    // The number of arguments that will be given to PETSc.
    int petscArgc;

    // The list of arguments that will be given to PETSc.
    char **petscArgv;

    // Value of the constant temperature in Kelvin
    double constTemp;

    // Name of the input temperature profile file.
    std::string tempProfileFileName;

    // Deal with the name of the network file.
    // @param arg Argument given for the network file.
    bool handleNetworkOption( std::string arg );

    // Callback when have seen the network option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handleNetworkOptionCB( Options* opts, std::string arg );

    // Deal with the constant temperature handler selection option.
    // @param arg Argument given to the constant temperature handler selection option.
    bool handleConstTemperatureOption( std::string arg);

    // Callback when have seen the constant temperature handlers option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handleConstTemperatureOptionCB( Options* opts, std::string arg);

    // Deal with the temperature file handler selection option.
    // @param arg Argument given to the temperature file handler selection option.
    bool handleTemperatureFileOption( std::string arg);

    // Callback when have seen the temperature file handlers option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handleTemperatureFileOptionCB( Options* opts, std::string arg);

    // Deal with the performance handler selection option.
    // @param arg Argument given to the handler selection option.
    bool handlePerfHandlersOption( std::string arg );

    // Callback when have seen the handlers option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handlePerfHandlersOptionCB( Options* opts, std::string arg );

    // Deal with the visualization handler selection option.
    // @param arg Argument given to the handler selection option.
    bool handleVizHandlersOption( std::string arg );

    // Callback when have seen the handlers option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handleVizHandlersOptionCB( Options* opts, std::string arg );

    // Deal with the PETSc arguments option.
    // @param arg Unused
    bool handlePetscOption( std::string arg );

    // Callback when have seen the PETSc arguments option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handlePetscOptionCB( Options* opts, std::string arg );

public:

    XolotlOptions( void );
    virtual ~XolotlOptions( void );

    // Read the parameters from the given file to set the different
    // xolotl options.
    // @param argc The number of arguments in the argv vector.
    // @param argv Vector of argument strings.
    virtual void readParams( int argc, char* argv[] );

    // Show our help message.
    // @param os The output stream upon which to print the help message.
    virtual void showHelp( std::ostream& os ) const;

    // Should we use const temperature handlers?
    // @return true if program should use const temp handlers
    bool useConstTemperatureHandlers( void ) const  { return useConstTempHandlers; }

    // Should we use temperature profile handlers?
    // @return true if program should use temperature profile handlers
    bool useTemperatureProfileHandlers( void ) const  { return useTempProfileHandlers; }

    // Obtain the name of the file containing the temperature profile data.
    // @return Name of the temperature file.
    std::string getTempProfileFilename( void ) const    { return tempProfileFileName; }

    // Obtain the value of the constant temperature to be used.
    // @return Constant temperature.
    double getConstTemperature( void ) const  { return constTemp; }

    // Which type of performance handlers should we use?
    // @return Enum value indicating type of performance handlers to use.
    xolotlPerf::IHandlerRegistry::RegistryType getPerfRegistryType(void) const { return perfRegistryType; }

    // Should we use the "standard" set of handlers?
    // If false, use dummy handlers.
    // @return true if program should use standard handlers, false if
    // should use dummy handlers.
    bool useVizStandardHandlers( void ) const  { return useVizStdHandlers; }

    // Obtain the name of the file holding the input network.
    // @return Name of the input network file.
    std::string getNetworkFilename( void ) const    { return networkFileName; }

    // Obtain the number of arguments for PETSc.
    // @return Number of arguments.
    int getPetscArgc( void ) const    { return petscArgc; }

    // Obtain the list of arguments for PETSc.
    // @return The list of arguments.
    char** getPetscArgv( void ) const    { return petscArgv; }
};

};

#endif // XOLOTLOPTIONS_H
