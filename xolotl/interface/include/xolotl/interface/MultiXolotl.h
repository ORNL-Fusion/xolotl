#pragma once

#include <memory>
#include <vector>

#include <xolotl/config.h>
#include <xolotl/interface/IXolotlInterface.h>
#include <xolotl/options/IOptions.h>
#include <xolotl/util/TimeStepper.h>

namespace xolotl
{
namespace interface
{
class ComputeContext;
class PetscContext;
class XolotlInterface;

class MultiXolotl : public IXolotlInterface
{
public:
	MultiXolotl(const std::shared_ptr<ComputeContext>& context,
		const std::shared_ptr<options::IOptions>& opts);

	virtual ~MultiXolotl();

	double
	previousTime() const noexcept;

	double
	currentTime() const noexcept;

	double
	currentDt() const noexcept;

	std::size_t
	currentStep() const noexcept;

	void
	solveXolotl() override;

	void
	solveStep();

private:
	void
	writeStopData();

	void
	startTimeStepper();

private:
	std::shared_ptr<ComputeContext> _computeContext;
	std::shared_ptr<options::IOptions> _options;
	util::TimeStepper _timeStepper;
	std::unique_ptr<PetscContext> _petscContext;
	std::unique_ptr<XolotlInterface> _primaryInstance;
	std::vector<std::unique_ptr<XolotlInterface>> _subInstances;
	std::vector<IdType> _subDOFs;
	std::vector<std::shared_ptr<RatesCapsule>> _constantRates;
	bool _restarting{false};
	bool _checkpointing{false};
	std::vector<std::string> _checkpointFiles;
};
} // namespace interface
} // namespace xolotl
