#ifndef OPTIONS_H
#define OPTIONS_H

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <memory>


namespace xolotlCore {

class Options
{
protected:
    // Information about a specific option we support
    struct OptInfo
    {
        // help message describing the option
        std::string helpMessage;

        // handler to call when option is seen
        // optArg might be an empty string if no argument is given
        bool (*optHandler)( Options* opts, std::string optArg );

        OptInfo(std::string _helpMessage,
                bool (*_handler)( Options* opts, std::string optArg ) )
          : helpMessage( _helpMessage ),
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

public:
    Options( void );
    virtual ~Options( void );


    // Read the parameters from the given file to set the different
    // xolotl options.
    //
    // @param argc The number of arguments in the argv vector.
    // @param argv Vector of argument strings.
    virtual void readParams( int argc, char* argv[] );


    // Read the lines in the parameter file.
    //
    // @return The vector of key and options
    virtual std::vector<std::string> ParamReadLine(std::shared_ptr<std::istream> inputstream);


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
