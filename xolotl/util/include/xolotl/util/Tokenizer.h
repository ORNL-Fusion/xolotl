#pragma once

#include <algorithm>
#include <sstream>
#include <vector>

#include <boost/tokenizer.hpp>

namespace xolotl
{
namespace util
{
namespace detail
{
template <typename T>
struct Extractor
{
	T
	operator()(const std::string& str)
	{
		T data;
		std::istringstream(str) >> data;
		return data;
	}
};

template <>
struct Extractor<std::string>
{
	const std::string&
	operator()(const std::string& str)
	{
		return str;
	}
};
} // namespace detail

template <typename T = std::string>
class Tokenizer
{
public:
	Tokenizer() = delete;

	explicit Tokenizer(const std::string& str, const char* delimiter = " ") :
		_tokens(str, boost::char_separator<char>(delimiter))
	{
	}

	void
	operator()(std::vector<T>& dataVector) const
	{
		std::transform(_tokens.begin(), _tokens.end(),
			std::back_inserter(dataVector), detail::Extractor<T>{});
	}

	std::vector<T>
	operator()() const
	{
		std::vector<T> vec;
		(*this)(vec);
		return vec;
	}

private:
	const boost::tokenizer<boost::char_separator<char>> _tokens;
};
} // namespace util
} // namespace xolotl
