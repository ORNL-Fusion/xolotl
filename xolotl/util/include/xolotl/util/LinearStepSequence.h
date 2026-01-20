#pragma once

#include <xolotl/util/StepSequence.h>

namespace xolotl
{
namespace util
{
class LinearStepSequence : public StepSequence
{
public:
	LinearStepSequence();

	LinearStepSequence(
		double initialValue, double finalValue, std::size_t rampSteps);

	virtual ~LinearStepSequence();

	void
	start() override;

	void
	step() override;

private:
	double _initialValue{};
	double _diff{};
	double _rampSteps{};
};
} // namespace util
} // namespace xolotl
