#pragma once

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
	setup() override;

	PetscErrorCode
	monitorLargest(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	startStop(TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	computeXenonRetention(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	computeHeliumRetention(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	computeAlloy(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	monitorScatter(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	monitorBubble(TS ts, PetscInt timestep, PetscReal time, Vec solution);

protected:
	std::shared_ptr<viz::IPlot> _scatterPlot;
};
} // namespace monitor
} // namespace solver
} // namespace xolotl
