#pragma once

#include <stdexcept>

namespace xolotl
{
namespace options
{
class InvalidOptionValue : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};
}
}
