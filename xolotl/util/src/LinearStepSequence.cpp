#include <xolotl/util/LinearStepSequence.h>

namespace xolotl
{
namespace util
{
LinearStepSequence::LinearStepSequence() = default;

LinearStepSequence::LinearStepSequence(
	double initialValue, double finalValue, std::size_t rampSteps) :
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
	auto newVal = _initialValue +
		_diff * (static_cast<double>(_currentStep) / _rampSteps);
	_currentValue = newVal;
}
} // namespace util
} // namespace xolotl
