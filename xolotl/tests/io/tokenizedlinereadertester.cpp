#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <tokenizedlinereader.h>

using namespace std;
using namespace xolotlCore;

BOOST_AUTO_TEST_SUITE(TokenizedLineReader_testSuite)

/**This operation checks default parsing setup of the TokenizedLineReader.*/
BOOST_AUTO_TEST_CASE(checkDefaultParsing) {

	// Local Declarations
	string doubleString = "0.0 1.0 5.0\r\n0.11 0.55 22.86 99.283\r\n0.000382 883.33 74.832\r\n";
	TokenizedLineReader<double> doubleReader;
	TokenizedLineReader<int> intReader;
	TokenizedLineReader<bool> boolReader;
	vector<double> dLine;

	//----- Check doubles -----//

	// Create the input stream
	shared_ptr<stringstream> testStream(new stringstream(stringstream::in | stringstream::out));
	*(testStream.get()) << doubleString;
	// Load the double reader
	doubleReader.setInputStream(testStream);
	// Get the first line and check it
	dLine = doubleReader.loadLine();
	BOOST_REQUIRE(!dLine.empty());
	BOOST_REQUIRE_EQUAL(3,dLine.size());
	BOOST_REQUIRE_CLOSE_FRACTION(0.0,dLine.at(0),0.001);
	BOOST_REQUIRE_CLOSE_FRACTION(1.0,dLine.at(1),0.001);
	BOOST_REQUIRE_CLOSE_FRACTION(5.0,dLine.at(2),0.001);
	// Get the second line and check it
	dLine = doubleReader.loadLine();
	BOOST_REQUIRE(!dLine.empty());
	BOOST_REQUIRE_EQUAL(4,dLine.size());
	BOOST_REQUIRE_CLOSE_FRACTION(0.11,dLine.at(0),0.001);
	BOOST_REQUIRE_CLOSE_FRACTION(0.55,dLine.at(1),0.001);
	BOOST_REQUIRE_CLOSE_FRACTION(22.86,dLine.at(2),0.001);
	BOOST_REQUIRE_CLOSE_FRACTION(99.283,dLine.at(3),0.0001);
	// Get the third line and check it
	dLine = doubleReader.loadLine();
	BOOST_REQUIRE(!dLine.empty());
	BOOST_REQUIRE_EQUAL(3,dLine.size());
	BOOST_REQUIRE_CLOSE_FRACTION(0.000382,dLine.at(0),0.001);
	BOOST_REQUIRE_CLOSE_FRACTION(883.33,dLine.at(1),0.001);
	BOOST_REQUIRE_CLOSE_FRACTION(74.382,dLine.at(2),0.0001);

	return;
}

/**This operation checks the TokenizedLineReader when the token is changed to a comma.*/
BOOST_AUTO_TEST_CASE(checkCSVParsing) {

	BOOST_FAIL( "Test not yet implemented" );

	return;
}

/**This operation checks the TokenizedLineReader when the comment token is changed.*/
BOOST_AUTO_TEST_CASE(checkCommentTokens) {

	BOOST_FAIL( "Test not yet implemented" );

	return;
}

BOOST_AUTO_TEST_SUITE_END()
