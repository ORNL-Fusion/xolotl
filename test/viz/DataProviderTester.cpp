#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <xolotl/viz/dataprovider/DataProvider.h>

using namespace std;
using namespace xolotl::viz;
using namespace dataprovider;

/**
 * This suite is responsible for testing the DataProvider.
 */
BOOST_AUTO_TEST_SUITE(DataProvider_testSuite)

/**
 * Method checking the ability to use a name.
 */
BOOST_AUTO_TEST_CASE(checkName) {

	// Create myDataProvider
	auto myDataProvider = make_shared<DataProvider>("myDataProvider");

	BOOST_REQUIRE_EQUAL("myDataProvider", myDataProvider->getName());
}

/**
 * Method checking you can add points to the data, get the data, and getDataMean().
 */
BOOST_AUTO_TEST_CASE(checkData) {

	// Create myDataProvider
	auto myDataProvider = make_shared<DataProvider>("myDataProvider");

	// Create a Point vector
	auto myPoints = make_shared< vector <Point> >();

	// And fill it with some Point
	Point aPoint;
	aPoint.value = 3.0; aPoint.t = 1.0; aPoint.x = 2.0;
	myPoints->push_back(aPoint);
	aPoint.value = 2.0; aPoint.t = 3.0; aPoint.x = 2.0;
	myPoints->push_back(aPoint);
	aPoint.value = 5.0; aPoint.t = 6.0; aPoint.x = -2.0;
	myPoints->push_back(aPoint);
	aPoint.value = -8.0; aPoint.t = 8.0; aPoint.x = 5.0;
	myPoints->push_back(aPoint);
	aPoint.value = 7.0; aPoint.t = 7.0; aPoint.x = 7.0;
	myPoints->push_back(aPoint);

	// Set these points in the myDataProvider
	myDataProvider->setPoints(myPoints);

	// Get them back
	auto dataPoints = myDataProvider->getDataPoints();

	// First check the size of the vector
	BOOST_REQUIRE_EQUAL(dataPoints->size(), myPoints->size());

	// Loop on all the points in dataPoints
	for (unsigned int i = 0; i < dataPoints->size(); i++) {

		// Check that all the fields are the same
		BOOST_REQUIRE_EQUAL(dataPoints->at(i).value, myPoints->at(i).value);
		BOOST_REQUIRE_EQUAL(dataPoints->at(i).t, myPoints->at(i).t);
		BOOST_REQUIRE_EQUAL(dataPoints->at(i).x, myPoints->at(i).x);
		BOOST_REQUIRE_EQUAL(dataPoints->at(i).y, myPoints->at(i).y);
		BOOST_REQUIRE_EQUAL(dataPoints->at(i).z, myPoints->at(i).z);
	}

	// Get the mean value of the data
	double mean = myDataProvider->getDataMean();

	// Check it is the right mean value:
	// (3 + 2 + 5 - 8 + 7) / 5
	BOOST_REQUIRE_EQUAL(mean, 1.8);
}

BOOST_AUTO_TEST_SUITE_END()
