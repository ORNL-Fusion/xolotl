#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <XolotlOptions.h>
#include <XolotlConfig.h>

using namespace std;
using namespace xolotlCore;

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
    xopts.readParams( fargc, fargv );

    // The XolotlOptions class has required positional parameters.
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
    xopts.readParams( fargc, fargv );

    // An unrecognized option should result in indicating
    // the program shouldn't run, and an error exit code.
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
    xopts.readParams( fargc, fargv );

    // If the help option is given,
    // the program shouldn't run
    // but it isn't a failure (so a success exit code).
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_SUCCESS );
}

BOOST_AUTO_TEST_CASE(badParamFileName)
{
    xolotlCore::XolotlOptions xopts;

	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/bla.txt");
	string filename = sourceDir + pathToFile;
    const char* fname = filename.c_str();

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
    xopts.readParams( fargc, fargv );

    // If the parameter file does not exist, xolotl should not run
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_FAILURE );
}

BOOST_AUTO_TEST_CASE(badParamFile)
{
    xolotlCore::XolotlOptions xopts;

	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/param_bad.txt");
	string filename = sourceDir + pathToFile;
    const char* fname = filename.c_str();

    // cons a command line with a bad option.
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    xopts.readParams( fargc, fargv );

    // An unrecognized parameter should result in indicating
    // the program shouldn't run, and an error exit code.
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_FAILURE );
}

BOOST_AUTO_TEST_CASE(goodParamFile)
{
    xolotlCore::XolotlOptions xopts;

	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/param_good.txt");
	string filename = sourceDir + pathToFile;
    const char* fname = filename.c_str();

    // cons a command line with a bad option.
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    xopts.readParams( fargc, fargv );

    // Xolotl should run with good parameters
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), true );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_SUCCESS );

    // Check the network filename
    BOOST_REQUIRE_EQUAL( xopts.getNetworkFilename(), "tungsten.txt" );

    // Check the temperature
    BOOST_REQUIRE_EQUAL( xopts.useConstTemperatureHandlers(), true );

    // Check the performance handler
    BOOST_REQUIRE_EQUAL( xopts.usePerfStandardHandlers(), true );

    // Check the performance handler
    BOOST_REQUIRE_EQUAL( xopts.useVizStandardHandlers(), true );

    // Check the PETSc options
    BOOST_REQUIRE_EQUAL( xopts.getPetscArgc(), 20 );
}

BOOST_AUTO_TEST_CASE(wrongPerfHandler)
{
    xolotlCore::XolotlOptions xopts;

	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/param_wrong.txt");
	string filename = sourceDir + pathToFile;
    const char* fname = filename.c_str();

    // cons a command line with a bad option.
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    xopts.readParams( fargc, fargv );

    // Xolotl should not be able to run with a wrong performance handler parameter
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_FAILURE );
}

BOOST_AUTO_TEST_CASE(bothTempArgs)
{
    xolotlCore::XolotlOptions xopts;

	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/param_bad_temp.txt");
	string filename = sourceDir + pathToFile;
    const char* fname = filename.c_str();

    // cons a command line with a bad option.
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = NULL;
    char** fargv = args;

    // attempt to parse the command line
    fargc -= 1;
    fargv += 1;
    xopts.readParams( fargc, fargv );

    // Xolotl should not be able to run when both temperature arguments are specified
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_FAILURE );
}

BOOST_AUTO_TEST_SUITE_END()



