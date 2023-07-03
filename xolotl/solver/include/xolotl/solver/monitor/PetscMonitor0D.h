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
	startStop(TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

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
	computeAlphaZr(TS ts, PetscInt timestep, PetscReal time, Vec solution);

	PetscErrorCode
	monitorBubble(TS ts, PetscInt timestep, PetscReal time, Vec solution);

protected:
	std::shared_ptr<viz::IPlot> _scatterPlot;

	std::vector<IdType> _clusterOrder;
};
} // namespace monitor
} // namespace solver
} // namespace xolotl
