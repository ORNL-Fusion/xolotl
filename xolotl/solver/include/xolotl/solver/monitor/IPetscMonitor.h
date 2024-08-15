#pragma once

#include <petscsys.h>
#include <petscts.h>

#include <string>

#include <xolotl/core/network/IReactionNetwork.h>
#include <xolotl/solver/monitor/IMonitor.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
class IPetscMonitor : public IMonitor
{
public:
	virtual ~IPetscMonitor()
	{
	}

	virtual void
	writeNetwork(MPI_Comm comm, const std::string& targetFileName,
		const std::string& srcFileName = "") = 0;

	virtual PetscErrorCode
	monitorTime(TS ts, PetscInt timestep, PetscReal time, Vec solution) = 0;

	virtual PetscErrorCode
	computeFluence(TS ts, PetscInt timestep, PetscReal time, Vec solution) = 0;

	virtual PetscErrorCode
	monitorPerf(TS ts, PetscInt timestep, PetscReal time, Vec solution) = 0;

	virtual PetscErrorCode
	monitorLargest(TS ts, PetscInt timestep, PetscReal time, Vec solution) = 0;

	virtual PetscErrorCode
	startStop(TS ts, PetscInt timestep, PetscReal time, Vec solution) = 0;

	virtual PetscErrorCode
	computeHeliumRetention(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) = 0;

	virtual PetscErrorCode
	computeXenonRetention(
		TS ts, PetscInt timestep, PetscReal time, Vec solution) = 0;

	virtual PetscErrorCode
	computeAlloy(TS ts, PetscInt timestep, PetscReal time, Vec solution) = 0;

	virtual PetscErrorCode
	computeT91(TS ts, PetscInt timestep, PetscReal time, Vec solution) = 0;

	virtual PetscErrorCode
	computeAlphaZr(TS ts, PetscInt timestep, PetscReal time, Vec solution) = 0;

	virtual PetscErrorCode
	monitorScatter(TS ts, PetscInt timestep, PetscReal time, Vec solution) = 0;

	virtual PetscErrorCode
	eventFunction(TS ts, PetscReal time, Vec solution, PetscScalar* fvalue) = 0;

	virtual PetscErrorCode
	postEventFunction(TS ts, PetscInt nevents, PetscInt eventList[],
		PetscReal time, Vec solution, PetscBool) = 0;
};
} // namespace monitor
} // namespace solver
} // namespace xolotl
