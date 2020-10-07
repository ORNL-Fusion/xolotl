#include <xolotl/perf/PerfHandlerRegistry.h>

namespace xolotl
{
namespace perf
{
std::shared_ptr<IPerfHandler> PerfHandlerRegistry::perfHandler;

void
PerfHandlerRegistry::set(const std::shared_ptr<IPerfHandler>& handler)
{
	perfHandler = handler;
}

std::shared_ptr<IPerfHandler>
PerfHandlerRegistry::get()
{
	if (!perfHandler) {
		throw std::runtime_error(
			"Request for performance handler before it has been registered");
	}
	return perfHandler;
}

// OpenMPI 2.x and 3.0 define its MPI_Datatype constants using a
// C-style cast to void*.  Clang++ objects to using these with
// in-class initializers and constexpr in our classes like ITimer.
// So we have to define them elsewhere (i.e., here).
const MPI_Datatype ITimer::MPIValType = MPI_DOUBLE;
const MPI_Datatype IEventCounter::MPIValType = MPI_UNSIGNED_LONG;
const MPI_Datatype IHardwareCounter::MPIValType = MPI_LONG_LONG_INT;
} // namespace perf
} // namespace xolotl
