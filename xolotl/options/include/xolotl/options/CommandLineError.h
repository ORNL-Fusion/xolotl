#pragma once

#include <stdexcept>

namespace xolotl
{
namespace options
{
class CommandLineError : public std::runtime_error
{
public:
	using std::runtime_error::runtime_error;
};
} // namespace options
} // namespace xolotl
