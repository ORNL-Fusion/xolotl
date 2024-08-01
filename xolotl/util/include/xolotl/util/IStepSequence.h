#pragma once

#include <cstdint>

namespace xolotl
{
namespace util
{
class IStepSequence
{
public:
	virtual ~IStepSequence()
	{
	}

	virtual void
	start() = 0;

	virtual void
	step() = 0;

	virtual double
	current() const noexcept = 0;

	virtual std::size_t
	currentStep() const noexcept = 0;

	virtual double
	at(std::size_t step) const = 0;

	virtual double
	partialSumAt(std::size_t step) const = 0;
};
} // namespace util
} // namespace xolotl
