#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "xolotlCore/commandline/XolotlOptions.h"



/**
 * Test suite for the Options class.
 */
BOOST_AUTO_TEST_SUITE (XolotlOptions_testSuite)

BOOST_AUTO_TEST_CASE(noOptions)
{
    xolotlCore::XolotlOptions xopts;

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

    // The XolotlOptions class has required positional parameters.
    BOOST_REQUIRE_EQUAL( nOptsUsed, 0 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_FAILURE );
}

BOOST_AUTO_TEST_CASE(badOption)
{
    xolotlCore::XolotlOptions xopts;

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
    BOOST_REQUIRE_EQUAL( nOptsUsed, 0 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_FAILURE );
}

BOOST_AUTO_TEST_CASE(helpOption)
{
    xolotlCore::XolotlOptions xopts;

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

BOOST_AUTO_TEST_CASE(networkFilename)
{
    xolotlCore::XolotlOptions xopts;
    const char* fname = "/home/joeuser/xolotl/tungsten_tiny.txt";

    // cons a command line with a valid network file name
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    int nOptsUsed = xopts.parseCommandLine( fargc, fargv );

    BOOST_REQUIRE_EQUAL( nOptsUsed+1, 1 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), true );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_SUCCESS );
    BOOST_REQUIRE_EQUAL( xopts.getNetworkFilename(), fname );
}

BOOST_AUTO_TEST_CASE(netFilenameAndBadOption)
{
    xolotlCore::XolotlOptions xopts;
    const char* fname = "/home/joeuser/xolotl/tungsten_tiny.txt";

    // cons a command line with a bad option.
    int fargc = 3;
    char* args[4];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = const_cast<char*>("--bogus");
    args[3] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    int nOptsUsed = xopts.parseCommandLine( fargc, fargv );

    // An unrecognized option should result in indicating
    // the program shouldn't run, and an error exit code.
    BOOST_REQUIRE_EQUAL( nOptsUsed+1, 2 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_FAILURE );
}

BOOST_AUTO_TEST_CASE(handlersNoArg)
{
    xolotlCore::XolotlOptions xopts;
    const char* fname = "/home/joeuser/xolotl/tungsten_tiny.txt";

    // cons a command line with a bad option.
    int fargc = 3;
    char* args[4];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = const_cast<char*>("--handlers");
    args[3] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    int nOptsUsed = xopts.parseCommandLine( fargc, fargv );

    // An unrecognized option should result in indicating
    // the program shouldn't run, and an error exit code.
    BOOST_REQUIRE_EQUAL( nOptsUsed+1, 2 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_FAILURE );
}


BOOST_AUTO_TEST_CASE(handlersDummyArg)
{
    xolotlCore::XolotlOptions xopts;
    const char* fname = "/home/joeuser/xolotl/tungsten_tiny.txt";

    // cons a command line
    int fargc = 4;
    char* args[5];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = const_cast<char*>("--handlers");
    args[3] = const_cast<char*>("dummy");
    args[4] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    int nOptsUsed = xopts.parseCommandLine( fargc, fargv );

    BOOST_REQUIRE_EQUAL( nOptsUsed+1, 3 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), true );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_SUCCESS );
    BOOST_REQUIRE_EQUAL( xopts.useStandardHandlers(), false );
}


BOOST_AUTO_TEST_CASE(handlersStdArg)
{
    xolotlCore::XolotlOptions xopts;
    const char* fname = "/home/joeuser/xolotl/tungsten_tiny.txt";

    // cons a command line 
    int fargc = 4;
    char* args[5];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = const_cast<char*>("--handlers");
    args[3] = const_cast<char*>("std");
    args[4] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    int nOptsUsed = xopts.parseCommandLine( fargc, fargv );

    BOOST_REQUIRE_EQUAL( nOptsUsed+1, 3 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), true );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_SUCCESS );
    BOOST_REQUIRE_EQUAL( xopts.useStandardHandlers(), true );
}


BOOST_AUTO_TEST_CASE(handlersBadArg)
{
    xolotlCore::XolotlOptions xopts;
    const char* fname = "/home/joeuser/xolotl/tungsten_tiny.txt";

    // cons a command line with a bad option.
    int fargc = 4;
    char* args[5];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = const_cast<char*>("--handlers");
    args[3] = const_cast<char*>("bogus");
    args[4] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    int nOptsUsed = xopts.parseCommandLine( fargc, fargv );

    // An unrecognized option should result in indicating
    // the program shouldn't run, and an error exit code.
    BOOST_REQUIRE_EQUAL( nOptsUsed+1, 3 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_FAILURE );
}


BOOST_AUTO_TEST_CASE(petscArg)
{
    xolotlCore::XolotlOptions xopts;
    const char* fname = "/home/joeuser/xolotl/tungsten_tiny.txt";

    // cons a command line
    int fargc = 6;
    char* args[7];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = const_cast<char*>("--petsc");
    args[3] = const_cast<char*>("one");
    args[4] = const_cast<char*>("two");
    args[5] = const_cast<char*>("three");
    args[6] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    int nOptsUsed = xopts.parseCommandLine( fargc, fargv );

    BOOST_REQUIRE_EQUAL( nOptsUsed+1, 2 );
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), true );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_SUCCESS );
    BOOST_REQUIRE_EQUAL( xopts.useStandardHandlers(), true );
}

BOOST_AUTO_TEST_SUITE_END()



