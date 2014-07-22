#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <Options.h>



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
    xopts.readParams( fargc, fargv );

    // The Option class does not like empty command line
    // a parameter file is always needed
    BOOST_REQUIRE_EQUAL( xopts.shouldRun(), false );
    BOOST_REQUIRE_EQUAL( xopts.getExitCode(), EXIT_FAILURE );
}

BOOST_AUTO_TEST_SUITE_END()

