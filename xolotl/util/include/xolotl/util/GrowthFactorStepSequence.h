#pragma once

#include <xolotl/util/IStepSequence.h>

namespace xolotl
{
namespace util
{
class GrowthFactorStepSequence : public IStepSequence
{
public:
	GrowthFactorStepSequence(
		double initialValue, double finalValue, double growthFactor = 1.0);

	virtual ~GrowthFactorStepSequence();

	void
	start() override;

	void
	step() override;

	double
	current() const noexcept override
	{
		return _currentValue;
	}

private:
	double _initialValue{};
	double _finalValue{};
	double _growthFactor{};
    std::size_t _currentStep{};
	double _currentValue{};
};
} // namespace util
} // namespace xolotl
