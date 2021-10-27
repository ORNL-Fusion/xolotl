#pragma once

#include <utility>
#include <vector>

#include <xolotl/perf/ITimer.h>
#include <xolotl/solver/monitor/PetscMonitor.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
class PetscMonitor2D : public PetscMonitor
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
	eventFunction(
		TS ts, PetscReal time, Vec solution, PetscScalar* fvalue) override;

	PetscErrorCode
	postEventFunction(TS ts, PetscInt nevents, PetscInt eventList[],
		PetscReal time, Vec solution, PetscBool) override;

	PetscErrorCode
	monitorSurface(TS ts, PetscInt timestep, PetscReal time, Vec solution);

protected:
	std::vector<std::vector<double>> _nSurf;
	std::vector<std::vector<double>> _nBulk;
	std::vector<std::vector<double>> _previousSurfFlux;
	std::vector<std::vector<double>> _previousBulkFlux;

	std::vector<std::pair<PetscInt, PetscInt>> _depthPositions;

	std::shared_ptr<viz::IPlot> _surfacePlot;

	// Timers
	std::shared_ptr<perf::ITimer> _gbTimer;
};
} // namespace monitor
} // namespace solver
} // namespace xolotl
