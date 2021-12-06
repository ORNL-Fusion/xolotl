#pragma once

#include <petscsys.h>
#include <petscts.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
extern double timeStepThreshold;

extern PetscErrorCode
checkTimeStep(TS ts);

extern PetscErrorCode
monitorTime(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);

extern PetscErrorCode
computeFluence(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);

extern PetscErrorCode
monitorPerf(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);

extern PetscErrorCode
monitorLargest(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);

extern PetscErrorCode
startStop(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);

extern PetscErrorCode
computeHeliumRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);

extern PetscErrorCode
computeXenonRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);

extern PetscErrorCode
computeAlloy(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);

extern PetscErrorCode
monitorScatter(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx);

extern PetscErrorCode
eventFunction(
	TS ts, PetscReal time, Vec solution, PetscScalar* fvalue, void* ctx);

extern PetscErrorCode
postEventFunction(TS ts, PetscInt nevents, PetscInt eventList[], PetscReal time,
	Vec solution, PetscBool fwd, void* ctx);
} // namespace monitor
} // namespace solver
} // namespace xolotl
