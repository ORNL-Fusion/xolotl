#include <algorithm>
#include <cmath>

#include <xolotl/util/GrowthFactorStepSequence.h>

namespace xolotl
{
namespace util
{
GrowthFactorStepSequence::GrowthFactorStepSequence(double initialValue,
	double finalValue, double growthFactor, std::size_t maxSteps) :
	StepSequence(maxSteps),
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
	_currentStep = 0;
	_currentValue = _initialValue;
}

void
GrowthFactorStepSequence::step()
{
	auto newVal = _currentValue * _growthFactor;
	newVal = std::min(newVal, _finalValue);
	_stepSize = newVal - _currentValue;
	_currentValue = newVal;
	++_currentStep;
}
} // namespace util
} // namespace xolotl
