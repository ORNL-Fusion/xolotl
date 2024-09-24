#include <cassert>

#include <xolotl/util/TimeStepper.h>

namespace xolotl
{
namespace util
{
TimeStepper::TimeStepper() = default;

TimeStepper::TimeStepper(std::unique_ptr<IStepSequence>&& stepSequence,
	double startTime, double endTime, std::size_t maxSteps) :
	_seq(std::move(stepSequence)),
	_maxSteps(maxSteps),
	_startTime(startTime),
	_endTime(endTime)
{
}

TimeStepper::TimeStepper(TimeStepper&& other) = default;
// {
// 	_seq(std::move(stepSequence)),
// 	_maxSteps(maxSteps),
// 	_startTime(startTime),
// 	_endTime(endTime)
// }
TimeStepper&
TimeStepper::operator=(TimeStepper&& other) = default;

void
TimeStepper::start()
{
	_previousTime = _startTime;
	_currentTime = _startTime;
	_seq->start();
}

void
TimeStepper::step()
{
	auto newTime = _currentTime + _seq->current();
	_seq->step();
	_previousTime = _currentTime;
	_currentTime = newTime;

	// TODO: should we add an option to require the end to be exact?
	//  if (_currentTime == _endTime) {
	//  	_currentTime = newTime;
	//  }
	//  else {
	//  	_currentTime = std::min(newTime, _endTime);
	//  }
}

[[nodiscard]] bool
TimeStepper::valid() const noexcept
{
	return _previousTime < _endTime &&
		(_maxSteps == 0 || _seq->currentStep() <= _maxSteps);
}

[[nodiscard]] double
TimeStepper::timeAtStep(std::size_t step) const
{
	assert(step <= _maxSteps);
	return _startTime + _seq->partialSumAt(step);
}

[[nodiscard]] double
TimeStepper::timeStepSizeAtStep(std::size_t step) const
{
	assert(step <= _maxSteps);
	return _seq->at(step);
}
} // namespace util
} // namespace xolotl
