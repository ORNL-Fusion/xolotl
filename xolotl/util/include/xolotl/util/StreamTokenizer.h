#pragma once

#include <iostream>

#include <xolotl/util/Tokenizer.h>

namespace xolotl
{
namespace util
{
class StreamTokenizer
{
public:
	StreamTokenizer() = delete;

	explicit StreamTokenizer(std::istream& inputStream,
		const char* dataDelimiter = " ", const char* commentDelimiter = "#") :
		_inputStream{inputStream},
		_delimiter{dataDelimiter},
		_commentDelimiter{commentDelimiter}
	{
	}

	void
	setDataDelimiter(const char* delim)
	{
		_delimiter = delim;
	}

	void
	setCommentDelimiter(const char* delim)
	{
		_commentDelimiter = delim;
	}

	template <typename T>
	void
	loadLine(std::vector<T>& dataVector) const
	{
		std::string line;
		if (!std::getline(_inputStream, line)) {
			std::cout << "StreamTokenizer: Stream is empty or closed."
					  << std::endl;
			return;
		}

		if (line.empty()) {
			return;
		}

		// If this line is a comment, skip it by calling this operation again
		// TODO: handle preceding whitespace and end-of-line comments
		if (line.find(_commentDelimiter) == 0) {
			loadLine(dataVector);
			return;
		}

		Tokenizer<T>{line, _delimiter}(dataVector);
	}

	template <typename T>
	std::vector<T>
	loadLine() const
	{
		std::vector<T> dataVec;
		loadLine(dataVec);
		return dataVec;
	}

private:
	std::istream& _inputStream;
	const char* _delimiter{" "};
	const char* _commentDelimiter{"#"};
};
} // namespace util
} // namespace xolotl
