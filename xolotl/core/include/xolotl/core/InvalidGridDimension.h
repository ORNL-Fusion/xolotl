#pragma once

#include <stdexcept>

namespace xolotl
{
namespace core
{
class InvalidGridDimension : public std::runtime_error
{
	using std::runtime_error::runtime_error;
};
} // namespace core
} // namespace xolotl
