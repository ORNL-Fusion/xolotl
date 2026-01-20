#pragma once

#include <iostream>
#include <string>

namespace xolotl
{
namespace util
{
struct Indent
{
	explicit Indent(std::size_t len = 4) : length(len), indent(indentString())
	{
	}

	const std::string
	indentString()
	{
		return std::string(length, ' ');
	}

	void
	increase()
	{
		indent += indentString();
	}

	void
	decrease()
	{
		indent.erase(indent.end() - length, indent.end());
	}

	Indent&
	operator++()
	{
		increase();
		return *this;
	}

	Indent&
	operator--()
	{
		decrease();
		return *this;
	}

	std::size_t length;
	std::string indent;
};

inline std::ostream&
operator<<(std::ostream& os, Indent indent)
{
	os << indent.indent;
	return os;
}
} // namespace util
} // namespace xolotl
