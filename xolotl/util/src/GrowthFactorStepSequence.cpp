#include <algorithm>
#include <cmath>

#include <xolotl/util/GrowthFactorStepSequence.h>

namespace xolotl
{
namespace util
{
GrowthFactorStepSequence::GrowthFactorStepSequence() = default;

GrowthFactorStepSequence::GrowthFactorStepSequence(double initialValue,
	double finalValue, double growthFactor, std::size_t initialStep) :
	StepSequence(initialStep),
	_initialValue(initialValue),
	_finalValue(finalValue),
	_growthFactor(growthFactor)
{
}

GrowthFactorStepSequence::~GrowthFactorStepSequence()
{
}

void
GrowthFactorStepSequence::start()
{
	_currentStep = _initialStep;
	_currentValue = _initialValue;
}

void
GrowthFactorStepSequence::step()
{
	auto newVal = _currentValue * _growthFactor;
	newVal = std::min(newVal, _finalValue);
	_currentValue = newVal;
	++_currentStep;
}

double
GrowthFactorStepSequence::at(std::size_t step) const
{
	return std::min(_initialValue * std::pow(_growthFactor, step), _finalValue);
}

double
GrowthFactorStepSequence::partialSumAt(std::size_t step) const
{
	auto val = _initialValue;
	auto sum = val;
	for (std::size_t i = _initialStep; i < step; ++i) {
		val = std::min(val * _growthFactor, _finalValue);
		sum += val;
	}
	return sum;
}
} // namespace util
} // namespace xolotl
