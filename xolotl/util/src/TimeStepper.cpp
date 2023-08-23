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

void
TimeStepper::start()
{
	_currentTime = _startTime;
	_seq->start();
}

void
TimeStepper::step()
{
	auto newTime = _currentTime + _seq->current();
	_seq->step();
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
	return _currentTime < _endTime &&
		(_maxSteps == 0 || _seq->currentStep() <= _maxSteps);
}
} // namespace util
} // namespace xolotl
