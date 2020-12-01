#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>

#include <xolotl/viz/LabelProvider.h>
#include <xolotl/viz/PlottingStyle.h>
#include <xolotl/viz/dataprovider/DataProvider.h>
#include <xolotl/viz/standard/plot/Plot.h>
#include <xolotl/viz/standard/plot/ScatterPlot.h>

using namespace std;
using namespace xolotl::viz;
using namespace dataprovider;
using namespace standard::plot;

/**
 * This suite is responsible for testing the Plot class.
 */
BOOST_AUTO_TEST_SUITE(Plot_testSuite)

/**
 * Method checking the ability to choose a PlottingStyle.
 */
BOOST_AUTO_TEST_CASE(checkPlottingStyle)
{
	// Create myScatterPlot
	auto myScatterPlot = make_shared<ScatterPlot>();

	PlottingStyle thePlottingStyle = LINE;

	// Set the PlottingStyle of myScatterPlot
	myScatterPlot->setPlottingStyle(thePlottingStyle);

	// Check it is the right one
	BOOST_REQUIRE_EQUAL(myScatterPlot->getPlottingStyle(), thePlottingStyle);
}

/**
 * Method checking everything related to the data provider.
 */
BOOST_AUTO_TEST_CASE(checkDataProvider)
{
	// Create myScatterPlot
	auto myScatterPlot = make_shared<ScatterPlot>();

	// Create myDataProvider
	auto myDataProvider = make_shared<DataProvider>();

	// Create a DataPoint vector
	auto myPoints = make_shared<vector<DataPoint>>();

	// And fill it with some Point
	DataPoint aPoint;
	aPoint.value = 3.0;
	aPoint.t = 1.0;
	aPoint.x = 2.0;
	myPoints->push_back(aPoint);
	aPoint.value = 2.0;
	aPoint.t = 3.0;
	aPoint.x = 2.0;
	myPoints->push_back(aPoint);
	aPoint.value = 5.0;
	aPoint.t = 6.0;
	aPoint.x = -2.0;
	myPoints->push_back(aPoint);
	aPoint.value = -8.0;
	aPoint.t = 8.0;
	aPoint.x = 5.0;
	myPoints->push_back(aPoint);
	aPoint.value = -7.0;
	aPoint.t = 7.0;
	aPoint.x = 7.0;
	myPoints->push_back(aPoint);

	// Set these points in the myDataProvider
	myDataProvider->setDataPoints(myPoints);

	// Set myDataProvider in myScatterPlot
	myScatterPlot->setDataProvider(myDataProvider);

	// Get the Points from the DataProvider from the Plot
	auto dataPoints = myScatterPlot->getDataProvider()->getDataPoints();

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

BOOST_AUTO_TEST_SUITE_END()
