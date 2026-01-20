#pragma once

#include <xolotl/util/IStepSequence.h>

namespace xolotl
{
namespace util
{
class StepSequence : public IStepSequence
{
public:
	StepSequence() = default;

	explicit StepSequence(std::size_t initialStep) : _initialStep(initialStep)
	{
	}

	virtual ~StepSequence()
	{
	}

	double
	current() const noexcept override
	{
		return _currentValue;
	}

	std::size_t
	currentStep() const noexcept override
	{
		return _currentStep;
	}

protected:
	std::size_t _initialStep{0};
	std::size_t _currentStep{};
	double _currentValue{};
};
} // namespace util
} // namespace xolotl
