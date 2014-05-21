#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "xolotlCore/commandline/Options.h"



/**
 * Test suite for the Options class.
 */
BOOST_AUTO_TEST_SUITE (Options_testSuite)

BOOST_AUTO_TEST_CASE(noOptions)
{
    xolotlCore::Options xopts;

    // cons a fake, empty command line.
    // Here, "empty" is what the shell would give us if no arguments
    // were provided (i.e., it contains the executable name or path),
    // and we skipped that executable name before calling the
    // parsing method as specified in its comment.
    int fargc = 1;
    char* args[2];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    int nOptsUsed = xopts.parseCommandLine( fargc, fargv );

    // As far as the Options class is concerned, an empty command line is OK.
    BOOST_REQUIRE_EQUAL( nOptsUsed+1, 0 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), true );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_SUCCESS );
}

BOOST_AUTO_TEST_CASE(badOption)
{
    xolotlCore::Options xopts;

    // cons a command line with a bad option.
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>("--bogus");
    args[2] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    int nOptsUsed = xopts.parseCommandLine( fargc, fargv );

    // An unrecognized option should result in indicating
    // the program shouldn't run, and an error exit code.
    BOOST_REQUIRE_EQUAL( nOptsUsed+1, 1 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_FAILURE );
}

BOOST_AUTO_TEST_CASE(helpOption)
{
    xolotlCore::Options xopts;

    // cons a command line with the help option
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>("--help");
    args[2] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    int nOptsUsed = xopts.parseCommandLine( fargc, fargv );

    // If the help option is given, 
    // the program shouldn't run
    // but it isn't a failure (so a success exit code).
    BOOST_REQUIRE_EQUAL( nOptsUsed+1, 1 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_SUCCESS );
}

BOOST_AUTO_TEST_SUITE_END()



