#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <tokenizedlinereader.h>

using namespace std;
using namespace xolotlCore;

// Test data
const string doubleString =
		"0.0 1.0 5.0\n0.11 0.55 22.86 99.283\n# Comment\n0.000382 883.33 74.832\n";
const string doubleCSVDoubleSlashCommentString =
		"0.0,1.0,5.0\n0.11,0.55,22.86,99.283\n// Comment\n0.000382,883.33,74.832\n";
const string intString = "1 3 5 7 9\n# Comment\n 0 2 4 6 \n";
const string intCSVString = "1,3,5,7,9\n# Comment\n0,2,4,6\n";

BOOST_AUTO_TEST_SUITE(TokenizedLineReader_testSuite)

/**This operation checks default parsing setup of the TokenizedLineReader.*/
BOOST_AUTO_TEST_CASE(checkDefaultParsing) {

	// Local Declarations
	TokenizedLineReader<double> doubleReader;
	TokenizedLineReader<int> intReader;
	vector<double> dLine;
	vector<int> iLine;

	//----- Check doubles -----//

	// Create the input stream
	shared_ptr<stringstream> doubleTestStream(
			new stringstream(stringstream::in | stringstream::out));
	*doubleTestStream << doubleString;
	// Load the double reader
	doubleReader.setInputStream(doubleTestStream);
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
	// The third line should skipped because it is a comment so, get the fourth
	// line and check it
	dLine = doubleReader.loadLine();
	BOOST_REQUIRE(!dLine.empty());
	BOOST_REQUIRE_EQUAL(3,dLine.size());
	BOOST_REQUIRE_CLOSE_FRACTION(0.000382,dLine.at(0),0.001);
	BOOST_REQUIRE_CLOSE_FRACTION(883.33,dLine.at(1),0.001);
	BOOST_REQUIRE_CLOSE_FRACTION(74.832,dLine.at(2),0.0001);

	//----- Check integers -----//

	// Create the input stream
	shared_ptr<stringstream> intTestStream(
			new stringstream(stringstream::in | stringstream::out));
	*intTestStream << intString;
	// Load the int reader
	intReader.setInputStream(intTestStream);
	// Get the first line and check it
	iLine = intReader.loadLine();
	BOOST_REQUIRE(!iLine.empty());
	BOOST_REQUIRE_EQUAL(5,iLine.size());
	BOOST_REQUIRE_EQUAL(1,iLine.at(0));
	BOOST_REQUIRE_EQUAL(3,iLine.at(1));
	BOOST_REQUIRE_EQUAL(5,iLine.at(2));
	BOOST_REQUIRE_EQUAL(7,iLine.at(3));
	BOOST_REQUIRE_EQUAL(9,iLine.at(4));
	// The second line should skipped because it is a comment so, get the fourth
	// line and check it
	iLine = intReader.loadLine();
	BOOST_REQUIRE(!iLine.empty());
	BOOST_REQUIRE_EQUAL(4,iLine.size());
	BOOST_REQUIRE_EQUAL(0,iLine.at(0));
	BOOST_REQUIRE_EQUAL(2,iLine.at(1));
	BOOST_REQUIRE_EQUAL(4,iLine.at(2));
	BOOST_REQUIRE_EQUAL(6,iLine.at(3));

	return;
}

/**This operation checks the TokenizedLineReader when the token is changed to a comma.*/
BOOST_AUTO_TEST_CASE(checkCSVParsing) {

	// Local Declarations
	TokenizedLineReader<int> intReader;
	vector<int> iLine;

	//----- Check integers -----//

	// Create the input stream
	shared_ptr<stringstream> intTestStream(
			new stringstream(stringstream::in | stringstream::out));
	*intTestStream << intCSVString;
	// Configure the delimiter
	intReader.setDelimiter(",");
	// Load the int reader
	intReader.setInputStream(intTestStream);
	// Get the first line and check it
	iLine = intReader.loadLine();
	BOOST_REQUIRE(!iLine.empty());
	BOOST_REQUIRE_EQUAL(5,iLine.size());
	BOOST_REQUIRE_EQUAL(1,iLine.at(0));
	BOOST_REQUIRE_EQUAL(3,iLine.at(1));
	BOOST_REQUIRE_EQUAL(5,iLine.at(2));
	BOOST_REQUIRE_EQUAL(7,iLine.at(3));
	BOOST_REQUIRE_EQUAL(9,iLine.at(4));
	// The second line should skipped because it is a comment so, get the fourth
	// line and check it
	iLine = intReader.loadLine();
	BOOST_REQUIRE(!iLine.empty());
	BOOST_REQUIRE_EQUAL(4,iLine.size());
	BOOST_REQUIRE_EQUAL(0,iLine.at(0));
	BOOST_REQUIRE_EQUAL(2,iLine.at(1));
	BOOST_REQUIRE_EQUAL(4,iLine.at(2));
	BOOST_REQUIRE_EQUAL(6,iLine.at(3));

	return;
}

/**This operation checks the TokenizedLineReader when the comment token is changed.*/
BOOST_AUTO_TEST_CASE(checkCommentDelimiter) {

	// Local Declarations
	TokenizedLineReader<double> doubleReader;
	vector<double> dLine;

	//----- Check doubles -----//

	// Create the input stream
	shared_ptr<stringstream> doubleTestStream(
			new stringstream(stringstream::in | stringstream::out));
	*doubleTestStream << doubleCSVDoubleSlashCommentString;
	// Set the delimiter
	doubleReader.setDelimiter(",");
	// Set the comment delimiter
	doubleReader.setCommentDelimiter("//");
	// Load the double reader
	doubleReader.setInputStream(doubleTestStream);
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
	// The third line should skipped because it is a comment so, get the fourth
	// line and check it
	dLine = doubleReader.loadLine();
	BOOST_REQUIRE(!dLine.empty());
	BOOST_REQUIRE_EQUAL(3,dLine.size());
	BOOST_REQUIRE_CLOSE_FRACTION(0.000382,dLine.at(0),0.001);
	BOOST_REQUIRE_CLOSE_FRACTION(883.33,dLine.at(1),0.001);
	BOOST_REQUIRE_CLOSE_FRACTION(74.832,dLine.at(2),0.0001);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
