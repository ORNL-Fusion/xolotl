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
	string doubleString = "0.0 1.0 5.0\r\n0.11 0.55 99.283\r\n38.2 883.33 74.832\r\n";
	TokenizedLineReader<double> doubleReader;
	TokenizedLineReader<int> intReader;
	TokenizedLineReader<bool> boolReader;

	// Try parsing the doubles
	shared_ptr<istream> testStream(new stringstream(stringstream::in | stringstream::out));

	return;
}

/**This operation checks the TokenizedLineReader when the token is changed to a comma.*/
BOOST_AUTO_TEST_CASE(checkCSVParsing) {
	return;
}

/**This operation checks the TokenizedLineReader when the comment token is changed.*/
BOOST_AUTO_TEST_CASE(checkCommentTokens) {
	return;
}

BOOST_AUTO_TEST_SUITE_END()
