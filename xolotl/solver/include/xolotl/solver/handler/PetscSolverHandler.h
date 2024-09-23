#ifndef PETSCSOLVERHANDLER_H
#define PETSCSOLVERHANDLER_H

// Includes
#include <xolotl/core/Types.h>
#include <xolotl/perf/IEventCounter.h>
#include <xolotl/perf/ITimer.h>
#include <xolotl/solver/handler/SolverHandler.h>
#include <xolotl/util/Profiling.h>

namespace xolotl
{
namespace solver
{
namespace handler
{
/**
 * This class and its subclasses realize the SolverHandler interface to solve
 * the diffusion-reaction problem with the PETSc solvers from Argonne
 * National Laboratory.
 *
 * This class does NOT implement most of the methods that are needed by the
 * PetscSolver. Only subclasses of this class must be used by the PetscSolver.
 */
class PetscSolverHandler : public SolverHandler
{
protected:
	std::optional<util::ProfileRegion> _profRegion;

	//! Partial derivatives at one grid point.
	Kokkos::View<double*> vals;

	//! Map of connectivities
	SparseFillMap dfill;

	//! The offset at the surface
	IdType surfaceOffset;

	/**
	 * A vector for holding the partial derivatives for one cluster in the order
	 * that PETSc expects. It is sized in the createSolverContext() operation.
	 *
	 * The vector is used for every cluster and immediately reset to zero before
	 * being used for the next. This allows the acquisition of the partial
	 * derivatives to take up minimal memory and require no additional dynamic
	 * allocations.
	 */
	std::vector<double> reactingPartialsForCluster;

	//! Times and counters
	std::shared_ptr<perf::ITimer> fluxTimer;
	std::shared_ptr<perf::ITimer> partialDerivativeTimer;
	std::shared_ptr<perf::IEventCounter> fluxCounter;
	std::shared_ptr<perf::IEventCounter> partialDerivativeCounter;

	/**
	 * Convert a C++ sparse fill map representation to the one that
	 * PETSc's DMDASetBlockFillsSparse() expects.
	 *
	 * @param dof The network's degrees of freedom.
	 * @param fillMap The C++ sparse file map to convert.
	 * @return The information from the fill map, in the format that
	 *      PETSc's DMDASetBlockFillsSparse() expects.
	 */
	static std::vector<PetscInt>
	ConvertToPetscSparseFillMap(size_t dof,
		const core::network::IReactionNetwork::SparseFillMap& fillMap);

	static std::array<std::vector<PetscInt>, 2>
	convertToCoordinateListPair(std::size_t dof,
		const core::network::IReactionNetwork::SparseFillMap& fillMap);

	static std::vector<core::RowColPair>
	convertToRowColPairList(std::size_t dof,
		const core::network::IReactionNetwork::SparseFillMap& fillMap);

public:
	/**
	 * Default constructor, deleted because we need to construct with objects.
	 */
	PetscSolverHandler() = delete;

	/**
	 * Construct a PetscSolverHandler.
	 *
	 * @param _network The reaction network to use.
	 * @param _perfHandler The perf handler to use.
	 */
	PetscSolverHandler(NetworkType& _network, perf::IPerfHandler& _perfHandler,
		const options::IOptions& options);

	virtual ~PetscSolverHandler();

	/**
	 * Reset all jacobian values to zero
	 */
	void
	resetJacobianValues();

	/**
	 * Set the number of grid points we want to move by at the surface.
	 * \see ISolverHandler.h
	 */
	void
	setSurfaceOffset(int offset)
	{
		surfaceOffset = offset;
		return;
	}
};
// end class PetscSolverHandler

} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */
#endif
