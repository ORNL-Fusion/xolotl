#pragma once

#include <memory>

#include <xolotl/util/IStepSequence.h>

namespace xolotl
{
namespace util
{
class TimeStepper
{
public:
	TimeStepper();

	TimeStepper(std::unique_ptr<IStepSequence>&& stepSequence, double startTime,
		double endTime, std::size_t maxSteps);

	void
	start();

	void
	step();

	TimeStepper&
	operator++()
	{
		step();
		return *this;
	}

	[[nodiscard]] bool
	valid() const noexcept;

	operator bool() const noexcept
	{
		return valid();
	}

	[[nodiscard]] double
	currentTime() const noexcept
	{
		return _currentTime;
	}

	[[nodiscard]] double
	previousTime() const noexcept
	{
		return _previousTime;
	}

	[[nodiscard]] double
	currentTimeStepSize() const noexcept
	{
		return _seq->current();
	}

	[[nodiscard]] std::size_t
	currentStep() const noexcept
	{
		return _seq->currentStep();
	}

private:
	std::unique_ptr<IStepSequence> _seq;
	std::size_t _maxSteps{};
	double _startTime{};
	double _endTime{};
	double _currentTime{};
	double _previousTime{};
};
} // namespace util
} // namespace xolotl
