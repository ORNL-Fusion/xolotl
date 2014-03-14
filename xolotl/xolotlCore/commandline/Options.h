#ifndef OPTIONS_H
#define OPTIONS_H

#include <iostream>
#include <string>
#include <map>


namespace xolotlCore {

class Options
{
protected:
    // Information about a specific option we support
    struct OptInfo
    {
        // is an argument required?
        bool argRequired;

        // help message describing the option
        std::string helpMessage;

        // handler to call when option is seen
        // optArg might be an empty string if no argument is given
        bool (*optHandler)( Options* opts, std::string optArg );

        OptInfo( bool _argRequired,
                    std::string _helpMessage,
                    bool (*_handler)( Options* opts, std::string optArg ) )
          : argRequired( _argRequired ),
            helpMessage( _helpMessage ),
            optHandler( _handler )
        { }
    };


    // Map of options we support, keyed by option switch string
    // (including leading dashes)
    typedef std::map<std::string, OptInfo*> OptionsMap;
    OptionsMap optionsMap;


    // Whether calling program should continue running after
    // parsing the command line.
    bool shouldRunFlag;


    // The exit code the program should use if we say the
    // program shouldn't run.
    int exitCode;


    // Deal with the help option.
    bool handleHelpOption( std::string arg );


    // Callback when have seen the help option.
    static bool handleHelpOptionCB( Options* opts, std::string optArg );

public:
    Options( void );
    virtual ~Options( void );


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


    // Should the program run after parsing the command line?
    // (Or did parsing the command line suggest the program 
    // shouldn't/couldn't run?)
    bool shouldRun( void ) const        { return shouldRunFlag; }


    // If program shouldn't run, what should its exit code be?
    int getExitCode( void ) const       { return exitCode; }
};

};

#endif // OPTIONS_H
