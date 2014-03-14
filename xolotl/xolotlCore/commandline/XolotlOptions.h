#ifndef XOLOTLOPTIONS_H
#define XOLOTLOPTIONS_H

#include <string>
#include "xolotlCore/commandline/Options.h"

namespace xolotlCore {

class XolotlOptions : public Options
{
private:
    // Use the "standard" set of handlers?
    bool useStdHandlers;

    // Name of the input network file.
    std::string netFileName;

    // Deal with the handler selection option.
    // @param arg Argument given to the handler selection option.
    bool handleHandlersOption( std::string arg );

    // Callback when have seen the handlers option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handleHandlersOptionCB( Options* opts, std::string arg );


    // Deal with the PETSc argument delimiter option.
    // @param arg Unused
    bool handlePetscOption( std::string arg );

    // Callback when have seen the PETSc argument delimiter option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handlePetscOptionCB( Options* opts, std::string arg );

public:
    XolotlOptions( void );


    // Parse the given command line for user-configurable settings.
    // We assume that the executable file/path has been skipped before
    // calling this method.  (E.g., the program's main() function
    // called this with something like 
    //   xopts.parseCommandLine( argc - 1, argv + 1 );
    //
    // @param argc The number of arguments in the argv vector.
    // @param argv Vector of argument strings.
    // @return Number of command line arguments used.
    virtual int parseCommandLine( int argc, char* argv[] );


    // Show our help message.
    // @param os The output stream upon which to print the help message.
    virtual void showHelp( std::ostream& os ) const;


    // Should we use the "standard" set of handlers?
    // If false, use dummy (stub) handlers.
    // TODO will we ever have more than {dummy, standard} such that 
    // we will need an enum?
    // @return true if program should use standard handlers, false if 
    // should use dummy handlers.
    bool useStandardHandlers( void ) const  { return useStdHandlers; }


    // Obtain the name of the file holding the input network.
    // @return Name of the input network file.
    std::string getNetworkFilename( void ) const    { return netFileName; }
};

};

#endif // XOLOTLOPTIONS_H
