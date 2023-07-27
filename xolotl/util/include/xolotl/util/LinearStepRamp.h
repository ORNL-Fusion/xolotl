#pragma once

#include <cstdint>

namespace xolotl
{
namespace util
{
class LinearStepRamp
{
public:
	LinearStepRamp() = default;

	LinearStepRamp(double initialValue, double finalValue, std::size_t steps) :
		_initialValue{initialValue},
		_diff{finalValue - _initialValue},
		_steps{static_cast<double>(steps - 1)}
	{
		start();
	}

	void
	start(double initialValue, double finalValue, std::size_t steps)
	{
		_initialValue = initialValue;
		_diff = finalValue - _initialValue;
		_steps = static_cast<double>(steps - 1);

		start();
	}

	void
	start()
	{
		_currentStep = 0.0;
		_currentValue = _initialValue;
	}

	void
	step()
	{
		++_currentStep;
		_currentValue = _initialValue + _diff * (_currentStep / _steps);
	}

	LinearStepRamp&
	operator++() noexcept
	{
		step();
		return *this;
	}

	double
	currentValue() const noexcept
	{
		return _currentValue;
	}

private:
	double _initialValue{};
	double _diff{};
	double _steps{};
	double _currentStep{};
	double _currentValue{};
};
} // namespace util
} // namespace xolotl
