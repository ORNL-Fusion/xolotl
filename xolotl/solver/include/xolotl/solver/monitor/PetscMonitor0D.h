#pragma once

#include <xolotl/perf/ITimer.h>
#include <xolotl/solver/monitor/PetscMonitor.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
class PetscMonitor0D : public PetscMonitor
{
public:
	using PetscMonitor::PetscMonitor;

	void
	setup(int loop) override;

	void
	keepFlux(std::vector<std::vector<std::vector<double>>>& nSurf,
		std::vector<std::vector<std::vector<double>>>& nBulk,
		std::vector<std::vector<std::vector<double>>>& surfFlux,
		std::vector<std::vector<std::vector<double>>>& bulkFlux) override
	{
	}

	void
	setFlux(std::vector<std::vector<std::vector<double>>>& nSurf,
		std::vector<std::vector<std::vector<double>>>& nBulk,
		std::vector<std::vector<std::vector<double>>>& surfFlux,
		std::vector<std::vector<std::vector<double>>>& bulkFlux) override
	{
	}

	PetscErrorCode
	monitorLargest(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	startStopImpl(TS ts, PetscInt timestep, PetscReal time, Vec solution,
		io::XFile& checkpointFile, io::XFile::TimestepGroup* tsGroup,
		const std::vector<std::string>& speciesNames) override;

	PetscErrorCode
	computeXenonRetention(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	computeAlloy(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	monitorScatter(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	computeAlphaZr(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	monitorBubble(TS ts, PetscInt timestep, PetscReal time, Vec solution);

	PetscErrorCode
	eventFunction(
		TS ts, PetscReal time, Vec solution, PetscScalar* fvalue) override;

	PetscErrorCode
	postEventFunction(TS ts, PetscInt nevents, PetscInt eventList[],
		PetscReal time, Vec solution, PetscBool) override;

protected:
	std::shared_ptr<viz::IPlot> _scatterPlot;

	std::vector<IdType> _clusterOrder;

	std::shared_ptr<perf::ITimer> _eventFuncTimer;
	std::shared_ptr<perf::ITimer> _postEventFuncTimer;

	double _previousInterI = 0.0;
	// V size threshold (44000 V)
	AmountType vThreshold = 44000;

	std::vector<std::vector<std::vector<std::pair<IdType, AmountType>>>>
		_vProductMap;
};
} // namespace monitor
} // namespace solver
} // namespace xolotl
