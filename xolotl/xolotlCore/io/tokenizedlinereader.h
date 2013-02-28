#ifndef TOKENIZEDLINEREADER_H
#define TOKENIZEDLINEREADER_H
//Begin section for file tokenizedlinereader.h
//End section for file tokenizedlinereader.h

using namespace std;

/**
 * This is a simple class for parsing lines from an input stream and
 * splitting that line based on a token. A line in this case is defined as any
 * sequence of characters that terminates with an end-of-line character such
 * as \r, \n or \r\n. The default token is a single space and the token can be
 * changed by calling the setToken() operation. Lines that are meant to be
 * comments will be skipped and this includes lines that start with a "#" by
 * default. The comment token can be changed by calling setCommentToken().
 * "In-line comments" where the comment starts at some later point in the line
 * are not allowed and will be parsed as data.</p><p>This class is templated so
 * that conversion to primitive types can be handled automatically for the
 * client class. This class is only meant to be templated with primitive types
 * and it will fail for classes or interfaces.
 */
template<class dataType>
class TokenizedLineReader {
	//Begin section for TokenizedLineReader
	//TODO: Add attributes that you want preserved
	//End section for TokenizedLineReader

private:

	//!The delimiter or "token" that separates elements in the line.
	const string dataToken;

	//!The character that represents a line that is a comment and should not be parsed.
	const string commentToken;

	//!The inputstream from which data should be read.
	istream inputstream;

public:

	//!The constructor
	TokenizedLineReader() :
		dataToken(" "), commentToken("#") {
		inputstream = NULL;
	}

	/**
	 * This operation sets the string (token) that represents the delimiter separating the different data elements in the line.
	 *
	 */
	void setToken(const string token) {
		//TODO Auto-generated method stub
		return;
	}

	/**
	 * This operation sets the character (token) that represents the comment character.
	 *
	 */
	void setCommentToken(const string cToken) {
		//TODO Auto-generated method stub
		return;
	}

	/**
	 * This operation sets the input stream that should be parsed.
	 */
	void setInputStream(istream stream) {
		inputstream = stream;
		return;
	}

	/**
	 * This operation loads the next line from the stream and returns an array
	 * of its elements by splitting on the token or null if there are no more
	 * lines in the stream.
	 */
	vector<dataType> loadLine() {
		//TODO Auto-generated method stub
	}

}; //end class TokenizedLineReader


#endif
