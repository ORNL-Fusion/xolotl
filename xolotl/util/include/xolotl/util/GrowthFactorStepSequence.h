#pragma once

#include <xolotl/util/StepSequence.h>

namespace xolotl
{
namespace util
{
class GrowthFactorStepSequence : public StepSequence
{
public:
	GrowthFactorStepSequence();

	GrowthFactorStepSequence(
		double initialValue, double finalValue, double growthFactor = 1.0);

	virtual ~GrowthFactorStepSequence();

	void
	start() override;

	void
	step() override;

private:
	double _initialValue{};
	double _finalValue{};
	double _growthFactor{};
};
} // namespace util
} // namespace xolotl
