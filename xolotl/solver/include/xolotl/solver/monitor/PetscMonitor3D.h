#pragma once

#include <array>
#include <vector>

#include <xolotl/solver/monitor/PetscMonitor.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
class PetscMonitor3D : public PetscMonitor
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
		nSurf = _nSurf;
		nBulk = _nBulk;
		surfFlux = _previousSurfFlux;
		bulkFlux = _previousBulkFlux;
	}

	void
	setFlux(std::vector<std::vector<std::vector<double>>>& nSurf,
		std::vector<std::vector<std::vector<double>>>& nBulk,
		std::vector<std::vector<std::vector<double>>>& surfFlux,
		std::vector<std::vector<std::vector<double>>>& bulkFlux) override
	{
		_nSurf = nSurf;
		_nBulk = nBulk;
		_previousSurfFlux = surfFlux;
		_previousBulkFlux = bulkFlux;
	}

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
	monitorSurfaceXY(TS ts, PetscInt timestep, PetscReal time, Vec solution);

	PetscErrorCode
	monitorSurfaceXZ(TS ts, PetscInt timestep, PetscReal time, Vec solution);

protected:
	std::vector<std::vector<std::vector<double>>> _nSurf;
	std::vector<std::vector<std::vector<double>>> _nBulk;
	std::vector<std::vector<std::vector<double>>> _previousSurfFlux;
	std::vector<std::vector<std::vector<double>>> _previousBulkFlux;

	std::vector<std::array<PetscInt, 3>> _depthPositions;

	std::shared_ptr<viz::IPlot> _surfacePlotXY;
	std::shared_ptr<viz::IPlot> _surfacePlotXZ;
};
} // namespace monitor
} // namespace solver
} // namespace xolotl
