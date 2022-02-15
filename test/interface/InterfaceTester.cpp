#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/framework.hpp>
#include <boost/test/unit_test.hpp>

#include <xolotl/interface/Interface.h>
#include <xolotl/test/CommandLine.h>

using namespace std;
using namespace xolotl;
using namespace interface;

/**
 * Test suite for the interface class.
 */
BOOST_AUTO_TEST_SUITE(Interface_testSuite)

BOOST_AUTO_TEST_CASE(simple0D)
{
	// Create the parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "vizHandler=dummy" << std::endl
			  << "petscArgs=-fieldsplit_0_pc_type redundant "
				 "-ts_max_snes_failures 200 "
				 "-pc_fieldsplit_detect_coupling "
				 "-ts_adapt_dt_max 10 "
				 "-pc_type fieldsplit "
				 "-fieldsplit_1_pc_type sor "
				 "-ts_max_time 1000 "
				 "-ts_max_steps 1 "
				 "-ts_exact_final_time stepover"
			  << std::endl
			  << "tempParam=900" << std::endl
			  << "perfHandler=dummy" << std::endl
			  << "flux=4.0e5" << std::endl
			  << "material=Fuel" << std::endl
			  << "dimensions=0" << std::endl
			  << "process=reaction" << std::endl
			  << "netParam=10 0 0 0 0" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};

	// Create and run the solver
	auto interface = xolotl::interface::XolotlInterface{cl.argc, cl.argv};
	interface.solveXolotl();

	// Get data to check
	auto concVector = interface.getConcVector();
	BOOST_REQUIRE_EQUAL(concVector[0][0][0][0].first, 0);
	BOOST_REQUIRE_SMALL(concVector[0][0][0][0].second, 1.0e-4);
	BOOST_REQUIRE_EQUAL(concVector[0][0][0][10].first, 10);
	BOOST_REQUIRE_EQUAL(concVector[0][0][0][10].second, 900);

	std::remove(parameterFile.c_str());
}

BOOST_AUTO_TEST_SUITE_END()
