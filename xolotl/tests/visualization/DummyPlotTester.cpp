#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <DummyPlot.h>
#include <PlottingStyle.h>
#include <DataProvider.h>

using namespace std;
using namespace xolotlViz;

/**
 * This suite is responsible for testing the DummyPlot class.
 */
BOOST_AUTO_TEST_SUITE(DummyPlot_testSuite)

/**
 * Method checking the non-ability to use a name.
 */
BOOST_AUTO_TEST_CASE(checkName) {

	// Create myDummyPlot
	auto myDummyPlot = make_shared<DummyPlot>("myDummyPlot");

	BOOST_REQUIRE_EQUAL("unused", myDummyPlot->getName());
}

/**
 * Method checking the non-ability to choose a PlottingStyle.
 */
BOOST_AUTO_TEST_CASE(checkPlottingStyle) {

	// Create myDummyPlot
	auto myDummyPlot = make_shared<DummyPlot>("myDummyPlot");

	PlottingStyle thePlottingStyle = LINE;

	// Set the PlottingStyle of myDummyPlot
	myDummyPlot->setPlottingStyle(thePlottingStyle);

	// Check it is the right one
	BOOST_REQUIRE_EQUAL(myDummyPlot->getPlottingStyle(), PlottingStyle());
}

/**
 * Method checking everything related to the data provider.
 */
BOOST_AUTO_TEST_CASE(checkDataProvider) {

	// Create myDummyPlot
	auto myDummyPlot = make_shared<DummyPlot>("myDummyPlot");

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
	aPoint.value = -7.0; aPoint.t = 7.0; aPoint.x = 7.0;
	myPoints->push_back(aPoint);

	// Set these points in the myDataProvider
	myDataProvider->setPoints(myPoints);

	// Set myDataProvider in myDummyPlot
	myDummyPlot->setDataProvider(myDataProvider);

	// Get the Points from the DataProvider from the Plot
	auto dataPoints = myDummyPlot->getDataProvider()->getDataPoints();

	// Loop on all the points in dataPoints
	for (unsigned int i = 0; i < dataPoints->size(); i++) {

		// Check that all the fields are the same
		BOOST_REQUIRE_EQUAL(dataPoints->at(i).value, myPoints->at(i).value);
		BOOST_REQUIRE_EQUAL(dataPoints->at(i).t, myPoints->at(i).t);
		BOOST_REQUIRE_EQUAL(dataPoints->at(i).x, myPoints->at(i).x);
		BOOST_REQUIRE_EQUAL(dataPoints->at(i).y, myPoints->at(i).y);
		BOOST_REQUIRE_EQUAL(dataPoints->at(i).z, myPoints->at(i).z);
	}
}

///**
// * Method checking the writing of the file.
// */
//BOOST_AUTO_TEST_CASE(checkWrite) {
//    BOOST_FAIL("checkWrite not implement yet");
//}

BOOST_AUTO_TEST_SUITE_END()
