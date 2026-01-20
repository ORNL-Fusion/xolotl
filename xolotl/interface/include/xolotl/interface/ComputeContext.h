#include <memory>

#include <mpi.h>

#include <Kokkos_Core.hpp>

#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace interface
{
class ComputeContext
{
public:
	ComputeContext(int& argc, const char* argv[]) :
		_kokkosContext((!Kokkos::is_initialized() && !Kokkos::is_finalized()) ?
				std::make_unique<Kokkos::ScopeGuard>(
					argc, const_cast<char**>(argv)) :
				nullptr)
	{
		if (!initialized()) {
			util::mpiInit(argc, argv);
			_mpiInitializedHere = true;
		}
	}

	~ComputeContext()
	{
		if (_mpiInitializedHere) {
			if (!finalized()) {
				MPI_Finalize();
			}
		}
	}

	static bool
	initialized()
	{
		int flag;
		MPI_Initialized(&flag);
		return flag != 0;
	}

	static bool
	finalized()
	{
		int flag;
		MPI_Finalized(&flag);
		return flag != 0;
	}

private:
	bool _mpiInitializedHere{false};
	std::unique_ptr<Kokkos::ScopeGuard> _kokkosContext;
};
} // namespace interface
} // namespace xolotl
