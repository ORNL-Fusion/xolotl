#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/framework.hpp>
#include <boost/test/unit_test.hpp>

#include <xolotl/interface/XolotlInterface.h>
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
	// Create a file with reaction data.
	std::ofstream reactionFile("reaction.dat");
	reactionFile
		<< "0 0 1 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "0 1 0 -0.160789024638881 -12.7387253625273 4.2488E+11 4.23317"
		<< std::endl
		<< "0 2 0 2.87060634573234 -8.37649260139129 9.9337E+11 3.43830674"
		<< std::endl
		<< "1 0 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 1 0 3.14296517286615 -1.4324963007117 0 0" << std::endl
		<< "1 2 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 3 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 4 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 5 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 6 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 7 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 8 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "Reactions" << std::endl
		<< "0 0 1 0.00596" << std::endl
		<< "0 1 0 0.00596" << std::endl
		<< "0 2 0 0.00596" << std::endl;
	reactionFile.close();

	// Create the parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "vizHandler=dummy" << std::endl
			  << "petscArgs=-fieldsplit_0_pc_type redundant "
				 "-ts_max_snes_failures 200 "
				 "-pc_fieldsplit_detect_coupling "
				 "-ts_adapt_dt_max 10 "
				 "-pc_type fieldsplit "
				 "-fieldsplit_1_pc_type jacobi "
				 "-ts_max_time 1000 "
				 "-ts_max_steps 1 "
				 "-ts_dt 0.1 "
				 "-ts_exact_final_time stepover"
			  << std::endl
			  << "tempParam=900" << std::endl
			  << "perfHandler=dummy" << std::endl
			  << "flux=4.0e5" << std::endl
			  << "material=Fuel" << std::endl
			  << "dimensions=0" << std::endl
			  << "process=reaction" << std::endl
			  << "reactionFilePath=reaction.dat" << std::endl
			  << "netParam=5 0 0 5 1" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};

	// Create and run the solver
	auto interface = xolotl::interface::XolotlInterface{cl.argc, cl.argv};
	interface.solveXolotl();

	// Get data to check
	auto concVector = interface.getConcVector();
	BOOST_REQUIRE_EQUAL(concVector[0][0][0][0].first, 0);
	BOOST_REQUIRE_CLOSE(concVector[0][0][0][0].second, 220899.426, 0.01);
	BOOST_REQUIRE_EQUAL(concVector[0][0][0][10].first, 10);
	BOOST_REQUIRE_CLOSE(
		concVector[0][0][0][10].second, -1.1602146071246332e-23, 0.01);

	std::remove(parameterFile.c_str());
	std::remove("reaction.dat");
}

BOOST_AUTO_TEST_SUITE_END()
