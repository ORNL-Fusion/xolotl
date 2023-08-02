#pragma once

#include <xolotl/util/IStepSequence.h>

namespace xolotl
{
namespace util
{
class LinearStepSequence : public IStepSequence
{
public:
	LinearStepSequence() = default;

	LinearStepSequence(
		double initialValue, double finalValue, std::size_t steps);

	void
	start() override;

	void
	step() override;

	LinearStepSequence&
	operator++() noexcept
	{
		step();
		return *this;
	}

	double
	current() const noexcept override
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
