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
	setup(int loop) override;

	void
	keepFlux(std::vector<std::vector<std::vector<double>>>& nSurf,
		std::vector<std::vector<std::vector<double>>>& nBulk,
		std::vector<std::vector<std::vector<double>>>& surfFlux,
		std::vector<std::vector<std::vector<double>>>& bulkFlux) override
	{
		nSurf.push_back(_nSurf);
		nBulk.push_back(_nBulk);
		surfFlux.push_back(_previousSurfFlux);
		bulkFlux.push_back(_previousBulkFlux);
	}

	void
	setFlux(std::vector<std::vector<std::vector<double>>>& nSurf,
		std::vector<std::vector<std::vector<double>>>& nBulk,
		std::vector<std::vector<std::vector<double>>>& surfFlux,
		std::vector<std::vector<std::vector<double>>>& bulkFlux) override
	{
		_nSurf = nSurf[0];
		_nBulk = nBulk[0];
		_previousSurfFlux = surfFlux[0];
		_previousBulkFlux = bulkFlux[0];
	}

	PetscErrorCode
	monitorLargest(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	startStopImpl(TS ts, PetscInt timestep, PetscReal time, Vec solution,
		io::XFile& checkpointFile, io::XFile::TimestepGroup* tsGroup,
		const std::vector<std::string>& speciesNames) override;

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
