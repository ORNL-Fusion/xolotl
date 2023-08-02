#include <xolotl/util/LinearStepSequence.h>

namespace xolotl
{
namespace util
{
LinearStepSequence::LinearStepSequence(
	double initialValue, double finalValue, std::size_t steps) :
	_initialValue{initialValue},
	_diff{finalValue - _initialValue},
	_steps{static_cast<double>(steps - 1)}
{
	start();
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
	_currentValue = _initialValue + _diff * (_currentStep / _steps);
}
} // namespace util
} // namespace xolotl
