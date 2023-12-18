#ifndef PETSCSOLVER2DHANDLER_H
#define PETSCSOLVER2DHANDLER_H

// Includes
#include <xolotl/solver/handler/PetscSolverHandler.h>

namespace xolotl
{
namespace solver
{
namespace handler
{
/**
 * This class is a subclass of PetscSolverHandler and implement all the methods
 * needed to solve the DR equations in 2D using PETSc from Argonne National
 * Laboratory.
 */
class PetscSolver2DHandler : public PetscSolverHandler
{
private:
	//! The position of the surface
	std::vector<IdType> surfacePosition;

public:
	PetscSolver2DHandler() = delete;

	/**
	 * Construct a PetscSolver2DHandler.
	 *
	 * @param _network The reaction network to use.
	 */
	PetscSolver2DHandler(
		NetworkType& _network, const options::IOptions& options) :
		PetscSolverHandler(_network, options)
	{
	}

	//! The Destructor
	~PetscSolver2DHandler()
	{
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	createSolverContext(DM& da);

	/**
	 * \see ISolverHandler.h
	 */
	void
	initializeSolverContext(DM& da, Mat& J);

	/**
	 * \see ISolverHandler.h
	 */
	void
	initializeConcentration(DM& da, Vec& C, DM& oldDA, Vec& oldC);

	/**
	 * \see ISolverHandler.h
	 */
	void
	initGBLocation(DM& da, Vec& C);

	/**
	 * \see ISolverHandler.h
	 */
	std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
	getConcVector(DM& da, Vec& C);

	/**
	 * \see ISolverHandler.h
	 */
	void
	setConcVector(DM& da, Vec& C,
		std::vector<
			std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>&
			concVector);

	/**
	 * \see ISolverHandler.h
	 */
	void
	updateConcentration(TS& ts, Vec& localC, Vec& F, PetscReal ftime);

	/**
	 * \see ISolverHandler.h
	 */
	void
	computeJacobian(TS& ts, Vec& localC, Mat& J, PetscReal ftime);

	/**
	 * \see ISolverHandler.h
	 */
	IdType
	getSurfacePosition(IdType j = -1, IdType k = -1) const
	{
		return surfacePosition[j];
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setSurfacePosition(IdType pos, IdType j = -1, IdType k = -1)
	{
		surfacePosition[j] = pos;

		return;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setSurfaceOffset(int offset, int j = -1, int k = -1)
	{
		return;
	}
};
// end class PetscSolver2DHandler

} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */
#endif
