#pragma once

#include <vector>

#include <xolotl/perf/ITimer.h>
#include <xolotl/solver/monitor/PetscMonitor.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
class PetscMonitor1D : public PetscMonitor
{
public:
	using PetscMonitor::PetscMonitor;

	void
	setup() override;

	PetscErrorCode
	monitorLargest(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	startStop(TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	computeHeliumRetention(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

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
	eventFunction(
		TS ts, PetscReal time, Vec solution, PetscScalar* fvalue) override;

	PetscErrorCode
	postEventFunction(TS ts, PetscInt nevents, PetscInt eventList[],
		PetscReal time, Vec solution, PetscBool) override;

	PetscErrorCode
	computeAlphaZr(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	checkNegative(TS ts, PetscInt timestep, PetscReal time, Vec solution);

	PetscErrorCode
	computeTRIDYN(TS ts, PetscInt timestep, PetscReal time, Vec solution);

	PetscErrorCode
	profileTemperature(TS ts, PetscInt timestep, PetscReal time, Vec solution);

	PetscErrorCode
	monitorSeries(TS ts, PetscInt timestep, PetscReal time, Vec solution);

protected:
	std::vector<double> _nSurf;
	std::vector<double> _nBulk;
	std::vector<double> _previousSurfFlux;
	std::vector<double> _previousBulkFlux;

	double _negThreshold = 0.0;

	std::vector<PetscInt> _depthPositions;

	std::shared_ptr<viz::IPlot> _seriesPlot;
	std::shared_ptr<viz::IPlot> _scatterPlot;

	// Timers
	std::shared_ptr<perf::ITimer> _initTimer;
	std::shared_ptr<perf::ITimer> _checkNegativeTimer;
	std::shared_ptr<perf::ITimer> _tridynTimer;
	std::shared_ptr<perf::ITimer> _startStopTimer;
	std::shared_ptr<perf::ITimer> _heRetentionTimer;
	std::shared_ptr<perf::ITimer> _xeRetentionTimer;
	std::shared_ptr<perf::ITimer> _scatterTimer;
	std::shared_ptr<perf::ITimer> _seriesTimer;
	std::shared_ptr<perf::ITimer> _eventFuncTimer;
	std::shared_ptr<perf::ITimer> _postEventFuncTimer;
};
} // namespace monitor
} // namespace solver
} // namespace xolotl
