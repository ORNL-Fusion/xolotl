#ifndef PETSCSOLVERHANDLER_H
#define PETSCSOLVERHANDLER_H

// Includes
#include "SolverHandler.h"

namespace xolotlSolver {

#ifndef CHECK_PETSC_ERROR
#define CHECK_PETSC_ERROR
/**
 * This operation checks a PETSc error code and throws an exception with given error message.
 *
 * @param errorCode The PETSc error code.
 * @param errMsg The error message in the thrown exception.
 */
inline void checkPetscError(PetscErrorCode errorCode, const char *errorMsg) {
	if (PetscUnlikely(errorCode))
		throw std::string(errorMsg);
}
#endif

/**
 * This class and its subclasses realize the ISolverHandler interface to solve the
 * advection-diffusion-reaction problem with the PETSc solvers from Argonne
 * National Laboratory.
 *
 * This class does NOT implement most of the methods that are needed by the
 * PetscSolver. Only subclasses of this class must be used by the PetscSolver.
 */
class PetscSolverHandler: public SolverHandler {
protected:

	//! Partial derivatives for all reactions at one grid point.
	Kokkos::View<double*> expVals;

	//! Map of connectivities
	SparseFillMap dfill;

	/**
	 * The last temperature on the grid. It is a vector to keep the temperature at each
	 * grid point but we know the temperature changes with depth only.
	 */
	std::vector<double> lastTemperature;

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

	/**
	 * Number of valid partial derivatives for each reactant.
	 */
	std::vector<PetscInt> reactionSize;

	/**
	 * Starting index of items for each reactant within the reactionIndices
	 * and reactionVals vectors.  E.g., the values for reactant i
	 * are located at
	 *      reactionIndices[reactionStartingIdx[i]+0],
	 *      reactionIndices[reactionStartingIdx[i]+1]
	 *      ...
	 *      reactionIndices[reactionStartingIdx[i]+reactionSize[i]-1]
	 */
	std::vector<size_t> reactionStartingIdx;

	/**
	 * Indices for partial derivatives for all the reactions at one grid point.
	 */
	std::vector<PetscInt> reactionIndices;

	/**
	 * Partial derivatives for all reactions at one grid point.
	 */
	std::vector<PetscScalar> reactionVals;

	//! Times and counters
	std::shared_ptr<xolotlPerf::ITimer> fluxTimer;
	std::shared_ptr<xolotlPerf::ITimer> partialDerivativeTimer;
	std::shared_ptr<xolotlPerf::IEventCounter> fluxCounter;
	std::shared_ptr<xolotlPerf::IEventCounter> partialDerivativeCounter;

	/**
	 * Convert a C++ sparse fill map representation to the one that
	 * PETSc's DMDASetBlockFillsSparse() expects.
	 *
	 * @param dof The network's degrees of freedom.
	 * @param fillMap The C++ sparse file map to convert.
	 * @return The information from the fill map, in the format that
	 *      PETSc's DMDASetBlockFillsSparse() expects.
	 */
	static std::vector<PetscInt> ConvertToPetscSparseFillMap(size_t dof,
			const xolotlCore::IReactionNetwork::SparseFillMap &fillMap);

public:

	/**
	 * Default constructor, deleted because we need to construct with objects.
	 */
	PetscSolverHandler() = delete;

	/**
	 * Construct a PetscSolverHandler.
	 *
	 * @param _network The reaction network to use.
	 */
	PetscSolverHandler(xolotlCore::IReactionNetwork &_network,
			NetworkType &_expNetwork) :
			SolverHandler(_network, _expNetwork), fluxTimer(
					xolotlPerf::getHandlerRegistry()->getTimer("Flux")), partialDerivativeTimer(
					xolotlPerf::getHandlerRegistry()->getTimer(
							"Partial Derivatives")), fluxCounter(
					xolotlPerf::getHandlerRegistry()->getEventCounter("Flux")), partialDerivativeCounter(
					xolotlPerf::getHandlerRegistry()->getEventCounter(
							"Partial Derivatives")) {
	}

};
//end class PetscSolverHandler

} /* end namespace xolotlSolver */
#endif
