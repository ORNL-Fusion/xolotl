#include <xolotl/util/LinearStepSequence.h>

namespace xolotl
{
namespace util
{
LinearStepSequence::LinearStepSequence(double initialValue, double finalValue,
	std::size_t rampSteps, std::size_t maxSteps) :
	StepSequence(maxSteps),
	_initialValue{initialValue},
	_diff{finalValue - _initialValue},
	_rampSteps{static_cast<double>(rampSteps)}
{
}

LinearStepSequence::~LinearStepSequence()
{
}

void
LinearStepSequence::start()
{
	_currentStep = 0.0;
	_currentValue = _initialValue;
}

void
LinearStepSequence::step()
{
	++_currentStep;
	if (_currentStep > _rampSteps) {
		return;
	}
	auto newVal = _initialValue + _diff * (_currentStep / _rampSteps);
	_stepSize = newVal - _currentValue;
	_currentValue = newVal;
}
} // namespace util
} // namespace xolotl
