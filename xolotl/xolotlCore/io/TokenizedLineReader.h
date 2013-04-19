#ifndef TOKENIZEDLINEREADER_H
#define TOKENIZEDLINEREADER_H

// Includes
#include <memory>
#include <iostream>
#include <sstream>

namespace xolotlCore {

/**
 * This is a simple class for parsing lines from an input stream and
 * splitting that line to retrieve delimiters. A line in this case is defined as any
 * sequence of characters that terminates with an end-of-line character such
 * as \r, \n or \r\n. The default delimiter is a single space and the delimiter can be
 * changed by calling the setDelimiter() operation. Delimiters found at the beginning
 * and end of the lines are removed. Lines that are meant to be
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
	std::string dataDelimiter;

	//!The character that represents a line that is a comment and should not be parsed.
	std::string commentDelimiter;

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

		dataDelimiter = delimiter;

		return;
	}

	/**
	 * This operation sets the character (delimiter) that represents the comment character.
	 *
	 */
	void setCommentDelimiter(const std::string cdelimiter) {

		commentDelimiter = cdelimiter;

		return;
	}

	/**
	 * This operation sets the input stream that should be parsed.
	 */
	void setInputStream(const std::shared_ptr<std::istream> stream) {
		inputstream = stream;
		return;
	}

	/**
	 * This operation loads the next line from the stream and returns an array
	 * of its elements by splitting on the delimiter or a vector of size zero
	 * if there are no more lines in the stream.
	 */
	std::vector<dataType> loadLine() {

		// Local Declarations
		size_t lastDelimiterPos = 0, nextDelimiterPos = 0, finalDelimiterPos =
				std::string::npos;
		std::vector<dataType> dataVector;
		std::string line, subLine;
		std::istringstream dataStream;
		dataType data;

		// Check the inputstream before reading from it
		if (!inputstream->eof() && inputstream->good()) {
			// Get the line
			std::getline(*inputstream, line);
			// Split it if it is not empty and does not start with the comment
			// character
			if (!line.empty()) {
				// If this line is a comment, skip it by calling this operation
				// again - FIXME! - Update to compiler to support C++11 and use
				// std::string.front()!
				if (line.find(commentDelimiter) == 0)
					return loadLine();
				//				// Remove delimiters at the beginning of the string
				if (line.find(dataDelimiter) == 0)
					line = line.substr(1);
				// Remove delimiters at the end of the string
				if (line.find(dataDelimiter, line.size() - 2) == line.size()
						- 1)
					line = line.erase(line.size() - 1);
				// Find the first instance of the delimiter
				nextDelimiterPos = line.find(dataDelimiter);
				// Only split the line if it contains the delimiter
				if (nextDelimiterPos != finalDelimiterPos) {
					// Walk across each piece of data in the line stopping only
					// when the end of the string is reached.
					while (lastDelimiterPos != finalDelimiterPos) {
						// Create the substring that holds the single piece of
						// data between the delimiters
						subLine = line.substr(lastDelimiterPos,
								nextDelimiterPos - lastDelimiterPos);
						// Clear the stream, reset the substring, load the data
						dataStream.clear();
						dataStream.str(subLine);
						dataStream >> data;
						dataVector.push_back(data);
						// Switch the delimiter positions and find the next
						lastDelimiterPos = (nextDelimiterPos
								!= finalDelimiterPos) ? nextDelimiterPos + 1
								: nextDelimiterPos;
						nextDelimiterPos = line.find(dataDelimiter,
								lastDelimiterPos);
					}
				} else {
					// Otherwise just put the whole line in the array
					dataStream >> data;
					dataVector.push_back(data);
				}
			}
		} else
			return dataVector;

		return dataVector;
	}

}; //end class TokenizedLineReader

} // end namespace xolotCore.io

#endif
