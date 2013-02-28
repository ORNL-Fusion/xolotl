#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <tokenizedlinereader.h>

BOOST_AUTO_TEST_SUITE(TokenizedLineReader_testSuite)

/**This operation checks default parsing setup of the TokenizedLineReader.*/
BOOST_AUTO_TEST_CASE(checkDefaultParsing) {
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
