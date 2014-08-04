#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <Options.h>
#include <XolotlConfig.h>

using namespace std;

/**
 * Test suite for the Options class.
 */
BOOST_AUTO_TEST_SUITE (Options_testSuite)

BOOST_AUTO_TEST_CASE(noOptions)
{
    xolotlCore::Options opts;

    // Build a fake, empty command line.
    // Here, "empty" is what the shell would give us if no arguments
    // were provided (i.e., it contains the executable name or path),
    // and we skipped that executable name before calling the
    // parsing method as specified in its comment.
    int fargc = 1;
    char* args[2];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = NULL;
    char** fargv = args;

    // Attempt to read the parameters
    fargc -= 1;
    fargv += 1;
    opts.readParams( fargc, fargv );

    // The Option class does not like empty command line
    // a parameter file is always needed
    BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
    BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);
}


BOOST_AUTO_TEST_CASE(badParamFileName)
{
    xolotlCore::Options opts;

	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/bla.txt");
	string filename = sourceDir + pathToFile;
    const char* fname = filename.c_str();

    // Build a command line with a non existing parameter file
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = NULL;
    char** fargv = args;

    // Attempt to read the parameters
    fargc -= 1;
    fargv += 1;
    opts.readParams(fargc, fargv);

    // If the parameter file does not exist, xolotl should not run
    BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
    BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);
}

BOOST_AUTO_TEST_CASE(badParamFile)
{
    xolotlCore::Options opts;

	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/param_bad.txt");
	string filename = sourceDir + pathToFile;
    const char* fname = filename.c_str();

    // Build a command line with a parameter file containing bad options
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = NULL;
    char** fargv = args;

    // Attempt to read the parameter file
    fargc -= 1;
    fargv += 1;
    opts.readParams(fargc, fargv);

    // An unrecognized parameter should result in indicating
    // the program shouldn't run, and an error exit code.
    BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
    BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);
}

BOOST_AUTO_TEST_CASE(goodParamFile)
{
    xolotlCore::Options opts;

	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/param_good.txt");
	string filename = sourceDir + pathToFile;
    const char* fname = filename.c_str();

    // Build a command line with a parameter file containing good options
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = NULL;
    char** fargv = args;

    // Attempt to read the parameter file
    fargc -= 1;
    fargv += 1;
    opts.readParams( fargc, fargv );

    // Xolotl should run with good parameters
    BOOST_REQUIRE_EQUAL(opts.shouldRun(), true);
    BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_SUCCESS);

    // Check the network filename
    BOOST_REQUIRE_EQUAL(opts.getNetworkFilename(), "tungsten.txt");

    // Check the temperature
    BOOST_REQUIRE_EQUAL(opts.useConstTemperatureHandlers(), true);
    BOOST_REQUIRE_EQUAL(opts.getConstTemperature(), 900.0);

    // Check if the maxHeFluence option is used
    BOOST_REQUIRE_EQUAL(opts.useMaxHeliumFluence(), true);
    BOOST_REQUIRE_EQUAL(opts.getMaxHeliumFluence(), 10.0);

    // Check if the heFlux option is used
    BOOST_REQUIRE_EQUAL(opts.useHeliumFlux(), true);
    BOOST_REQUIRE_EQUAL(opts.getHeliumFlux(), 1.5);

    // Check the performance handler
    BOOST_REQUIRE_EQUAL(opts.usePerfStandardHandlers(), true);

    // Check the performance handler
    BOOST_REQUIRE_EQUAL(opts.useVizStandardHandlers(), true);

    // Check the PETSc options
    BOOST_REQUIRE_EQUAL(opts.getPetscArgc(), 20);
}

BOOST_AUTO_TEST_CASE(wrongPerfHandler)
{
    xolotlCore::Options opts;

	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/param_wrong.txt");
	string filename = sourceDir + pathToFile;
    const char* fname = filename.c_str();

    // Build a command line with a parameter file containing a wrong performance handler option
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = NULL;
    char** fargv = args;

    // Attempt to read the parameter file
    fargc -= 1;
    fargv += 1;
    opts.readParams(fargc, fargv);

    // Xolotl should not be able to run with a wrong performance handler parameter
    BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
    BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);
}

BOOST_AUTO_TEST_CASE(bothTempArgs)
{
    xolotlCore::Options opts;

	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/param_bad_temp.txt");
	string filename = sourceDir + pathToFile;
    const char* fname = filename.c_str();

    // Build a command line with a parameter file containing the two temperature options
    int fargc = 2;
    char* args[3];
    args[0] = const_cast<char*>("./xolotl");
    args[1] = const_cast<char*>(fname);
    args[2] = NULL;
    char** fargv = args;

    // Attempt to read the parameter file
    fargc -= 1;
    fargv += 1;
    opts.readParams(fargc, fargv);

    // Xolotl should not be able to run when both temperature arguments are specified
    BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
    BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);
}

BOOST_AUTO_TEST_SUITE_END()
