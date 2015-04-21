#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include <TemperatureProfileHandler.h>

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the TemperatureHandler.
 */
BOOST_AUTO_TEST_SUITE (TemperatureProfileHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(check_getTemperature) {
	// Create a file with temperature profile data
	// First column with the time and the second with
	// the temperature at that time.
	std::ofstream writeTempFile("tempFile.dat");
	writeTempFile << "0.0 2.0 \n"
			"1.0 1.99219766723 \n"
			"2.0 1.87758256189 \n"
			"3.0 1.4311765168 \n"
			"4.0 0.583853163453 \n"
			"5.0 0.000137654918313 \n"
			"6.0 0.789204200569 \n"
			"7.0 1.9875147713 \n"
			"8.0 0.854499966191 \n"
			"9.0 0.235300873168 \n"
			"10.0 1.99779827918";
	writeTempFile.close();

	// Create and initialize the temperature profile handler
	auto testTemp = make_shared<TemperatureProfileHandler>("tempFile.dat");
	testTemp->initializeTemperature();
	std::vector<double> pos = {1.142857142857143, 0.0, 0.0};

	// Vector to hold the user defined time values
	std::vector<double> t;
	t.push_back(1.57894736842);
	t.push_back(4.21052631579);
	t.push_back(7.89473684211);
	t.push_back(2.10526315789);
	t.push_back(-1.0);
	t.push_back(12.0);

	// Vector to hold the true interpolated values
	std::vector<double> trueInterp;
	trueInterp.push_back(1.925841553612);
	trueInterp.push_back(0.460965687972);
	trueInterp.push_back(0.973764682518);
	trueInterp.push_back(1.83059245188);
	trueInterp.push_back(2.0);
	trueInterp.push_back(1.997798279179);

	// Vector to hold interpolated values
	std::vector<double> tempInterp;
	for (int i = 0; i < t.size(); i++) {
		tempInterp.push_back(testTemp->getTemperature(pos, t[i]));
	}

	// Verify the values
	for(int j = 0; j < t.size(); j++)
		BOOST_REQUIRE_CLOSE(tempInterp[j], trueInterp[j], 10e-8);

    // Remove the created file
	std::string tempFile = "tempFile.dat";
    std::remove(tempFile.c_str());

    return;
}

BOOST_AUTO_TEST_SUITE_END()
