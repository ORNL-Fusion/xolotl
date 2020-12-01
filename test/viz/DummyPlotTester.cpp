#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>

#include <xolotl/viz/PlottingStyle.h>
#include <xolotl/viz/dataprovider/DataProvider.h>
#include <xolotl/viz/dummy/DummyPlot.h>

using namespace std;
using namespace xolotl::viz;
using namespace dataprovider;
using namespace dummy;

/**
 * This suite is responsible for testing the DummyPlot class.
 */
BOOST_AUTO_TEST_SUITE(DummyPlot_testSuite)

/**
 * Method checking the non-ability to choose a PlottingStyle.
 */
BOOST_AUTO_TEST_CASE(checkPlottingStyle)
{
	// Create myDummyPlot
	auto myDummyPlot = make_shared<DummyPlot>();

	PlottingStyle thePlottingStyle = LINE;

	// Set the PlottingStyle of myDummyPlot
	myDummyPlot->setPlottingStyle(thePlottingStyle);

	// Check it is the right one
	BOOST_REQUIRE_EQUAL(myDummyPlot->getPlottingStyle(), PlottingStyle());
}

BOOST_AUTO_TEST_SUITE_END()
