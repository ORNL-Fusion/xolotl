#ifndef TOKENIZEDLINEREADER_H
#define TOKENIZEDLINEREADER_H

namespace xolotlCore {

/**
 * This is a simple class for parsing lines from an input stream and
 * splitting that line to retrieve delimiters. A line in this case is defined as any
 * sequence of characters that terminates with an end-of-line character such
 * as \r, \n or \r\n. The default delimiter is a single space and the delimiter can be
 * changed by calling the setDelimiter() operation. Lines that are meant to be
 * comments will be skipped and this includes lines that start with a "#" by
 * default. The comment delimiter can be changed by calling setCommentDelimiter().
 * "In-line comments" where the comment starts at some later point in the line
 * are not allowed and will be parsed as data. This class is templated so
 * that conversion to primitive types can be handled automatically for the
 * client class. This class is only meant to be templated with primitive types
 * and it will fail for classes or interfaces.
 */
template<class dataType>
class TokenizedLineReader {

private:

	//!The delimiter or "delimiter" that separates elements in the line.
	const std::string dataDelimiter;

	//!The character that represents a line that is a comment and should not be parsed.
	const std::string commentDelimiter;

	//!The inputstream from which data should be read.
	std::shared_ptr<std::istream> inputstream;

public:

	//!The constructor
	TokenizedLineReader() :
		dataDelimiter(" "), commentDelimiter("#") {
	}

	/**
	 * This operation sets the string (delimiter) that represents the delimiter separating the different data elements in the line.
	 *
	 */
	void setDelimiter(const std::string delimiter) {
		//TODO Auto-generated method stub
		return;
	}

	/**
	 * This operation sets the character (delimiter) that represents the comment character.
	 *
	 */
	void setCommentDelimiter(const std::string cdelimiter) {
		//TODO Auto-generated method stub
		return;
	}

	/**
	 * This operation sets the input stream that should be parsed.
	 */
	void setInputStream(std::shared_ptr<std::istream> stream) {
		inputstream = stream;
		return;
	}

	/**
	 * This operation loads the next line from the stream and returns an array
	 * of its elements by splitting on the delimiter or null if there are no more
	 * lines in the stream.
	 */
	std::vector<dataType> loadLine() {

		// Local Declarations
		std::vector<dataType> lineVector;
		std::string line;
		dataType dataLine;

		// Check the inputstream before reading from it
		if (!inputstream->eof() && inputstream->good()) {
			// Get the line
			std::getline(*(inputstream.get()), line);
			// Split it if it is not empty and does not start with the comment character
			std::cout << line << std::endl;
			if (!line.empty() && line.begin() == commentDelimiter) {
				// Create a stream from which to read the string
				std::istringstream dataStream(line);
				// Only split the line if it contains the delimiter
				if (line.find(dataDelimiter) != std::string::npos) {
					dataStream >> dataLine;
					lineVector.push_back(dataLine);
				} else
					// Otherwise just put the whole line in the array
					dataStream >> dataLine;
					lineVector.push_back(dataLine);
			}
		}

		return lineVector;
	}

}; //end class TokenizedLineReader

} // end namespace xolotCore.io

#endif
