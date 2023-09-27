#pragma once

#include <cstdint>

namespace xolotl
{
namespace util
{
class IStepSequence
{
public:
	virtual void
	start() = 0;

	virtual void
	step() = 0;

	virtual double
	current() const noexcept = 0;

	virtual std::size_t
	currentStep() const noexcept = 0;
};
} // namespace util
} // namespace xolotl
