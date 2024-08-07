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

	GrowthFactorStepSequence(double initialValue, double finalValue,
		double growthFactor = 1.0, std::size_t initialStep = 0);

	virtual ~GrowthFactorStepSequence();

	void
	start() override;

	void
	step() override;

	double
	at(std::size_t step) const override;

	double
	partialSumAt(std::size_t step) const override;

private:
	double _initialValue{};
	double _finalValue{};
	double _growthFactor{};
};
} // namespace util
} // namespace xolotl
