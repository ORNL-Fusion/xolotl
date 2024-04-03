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

BOOST_AUTO_TEST_CASE(compareOptions)
{
	auto cases = std::vector<std::string>{"benchmark_AZr_1", "benchmark_NE_1",
		"benchmark_NE_2", "benchmark_NE_3", "benchmark_NE_4", "benchmark_NE_5",
		"benchmark_PSI_10", "benchmark_PSI_1", "benchmark_PSI_2",
		"benchmark_PSI_3", "benchmark_PSI_4", "benchmark_PSI_5",
		"benchmark_PSI_7", "benchmark_PSI_8", "benchmark_PSI_9", "multi_0D",
		"multi_1D", "multi", "single_1D", "system_Alloy_1", "system_AZr_0",
		"system_Fe_1", "system_NE_0", "system_NE_1", "system_NE_2",
		"system_NE_3", "system_NE_4", "system_PSI_1", "system_PSI_2",
		"system_PSI_3", "system_PSI_4", "system_PSI_5"};

	::xolotl::fs::current_path(dataDir);

	auto readAndPrint = [](std::shared_ptr<IOptions> opts,
							const std::string& fileName) {
		const char* argv[] = {"", fileName.c_str()};
		opts->readParams(2, argv);
		std::ostringstream oss;
		opts->printAll(oss);
		return oss.str();
	};

	for (auto&& caseName : cases) {
		auto fileBaseName = "params_" + caseName;

		auto confStr = readAndPrint(
			std::make_shared<ConfOptions>(), fileBaseName + ".txt");
		auto jsonStr = readAndPrint(
			std::make_shared<JSONOptions>(), fileBaseName + ".json");

		BOOST_REQUIRE(confStr == jsonStr);
	}
}

BOOST_AUTO_TEST_SUITE_END()
