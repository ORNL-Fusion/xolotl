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

	StepSequence(std::size_t maxSteps) : _maxSteps(maxSteps)
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

	double
	stepSize() const noexcept override
	{
		return _stepSize;
	}

	bool
	valid() const noexcept override
	{
		return _maxSteps == 0 || _currentStep <= _maxSteps;
	}

protected:
	std::size_t _maxSteps{};
	double _currentStep{};
	double _currentValue{};
	double _stepSize{};
};
} // namespace util
} // namespace xolotl
