#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/program_options/errors.hpp>
#include <boost/test/framework.hpp>
#include <boost/test/unit_test.hpp>
namespace bpo = boost::program_options;

#include <xolotl/options/ConfOptions.h>
#include <xolotl/options/InvalidOptionValue.h>
#include <xolotl/options/JSONOptions.h>
#include <xolotl/options/Options.h>
#include <xolotl/test/config.h>
const std::string dataDir = TO_STRING(XOLOTL_TEST_DATA_DIR);

using namespace std;
using namespace xolotl::options;

/**
 * Test suite for the Options class.
 */
BOOST_AUTO_TEST_SUITE(Options_testSuite)

BOOST_AUTO_TEST_CASE(createOptions)
{
	const char* argv[] = {"", "test.txt"};
	auto opts = ::xolotl::options::createOptions(2, argv);
	BOOST_REQUIRE(dynamic_cast<ConfOptions*>(opts.get()));
	argv[1] = "test.json";
	opts = ::xolotl::options::createOptions(2, argv);
	BOOST_REQUIRE(dynamic_cast<JSONOptions*>(opts.get()));
	argv[1] = "test.conf";
	opts = ::xolotl::options::createOptions(2, argv);
	BOOST_REQUIRE(dynamic_cast<ConfOptions*>(opts.get()));
}

BOOST_AUTO_TEST_SUITE_END()
