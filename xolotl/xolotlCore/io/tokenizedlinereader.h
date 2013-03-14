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
			std::getline(*(inputstream.get()), line);
			// Split it if it is not empty and does not start with the comment character
			std::cout << line << std::endl;
			if (!line.empty()) {
				// If this line is a comment, skip it by calling this operation again
				// FIXME! - Update to compiler to support C++11 and use std::string.front()!
				if (line.find(commentDelimiter) == 0)
					return loadLine();
				// Only split the line if it contains the delimiter
				nextDelimiterPos = line.find(dataDelimiter);
				if (nextDelimiterPos != finalDelimiterPos) {
					while (lastDelimiterPos != finalDelimiterPos) {
						std::cout << "last and next delimiter positions "
								<< lastDelimiterPos + 1 << " "
								<< nextDelimiterPos << std::endl;
						subLine = line.substr(lastDelimiterPos + 1,
								nextDelimiterPos - lastDelimiterPos + 1);
						// SUBLINE is not updating correctly
						dataStream.clear();
						dataStream.str(subLine);
						// AND dataLine is wrong
						dataStream >> data;
						std::cout << "subLine = " << subLine << " with data "
								<< data << std::endl;
						dataVector.push_back(data);
						lastDelimiterPos = nextDelimiterPos;
						nextDelimiterPos = line.find(dataDelimiter,
								lastDelimiterPos + 1);
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
