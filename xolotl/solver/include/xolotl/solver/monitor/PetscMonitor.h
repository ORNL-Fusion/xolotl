#pragma once

#include <memory>

#include <xolotl/solver/handler/ISolverHandler.h>
#include <xolotl/solver/monitor/IPetscMonitor.h>
#include <xolotl/viz/IPlot.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
class PetscMonitor : public IPetscMonitor
{
public:
	PetscMonitor();

	PetscMonitor(
		TS ts, const std::shared_ptr<handler::ISolverHandler>& solverHandler);

	virtual ~PetscMonitor();

	void
	writeNetwork(MPI_Comm comm, const std::string& targetFileName,
		const std::string& srcFileName = "") override;

	PetscErrorCode
	monitorTime(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	computeFluence(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	monitorPerf(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	computeHeliumRetention(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) override;

	PetscErrorCode
	computeAlphaZr(
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

protected:
	TS _ts;

	std::shared_ptr<handler::ISolverHandler> _solverHandler;

	std::shared_ptr<viz::IPlot> _perfPlot;

	std::vector<IdType> _iClusterIds;
	int _largestClusterId = -1;
	double _largestThreshold = 1.0e-12;

	PetscInt _previousTSNumber = -1;

	double _sputteringYield = 0.0;

	double _nHeliumBurst = 0.0;
	double _nDeuteriumBurst = 0.0;
	double _nTritiumBurst = 0.0;

	PetscReal _hdf5Stride = 0.0;
	PetscInt _hdf5Previous = 0;
	std::string _hdf5OutputName = "xolotlStop.h5";
};
} // namespace monitor
} // namespace solver
} // namespace xolotl
