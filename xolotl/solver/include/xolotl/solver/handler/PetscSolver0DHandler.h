#ifndef PETSCSOLVER0DHANDLER_H
#define PETSCSOLVER0DHANDLER_H

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
 * needed to solve the DR equations in 0D using PETSc from Argonne National
 * Laboratory.
 */
class PetscSolver0DHandler : public PetscSolverHandler
{
public:
	PetscSolver0DHandler() = delete;

	/**
	 * Construct a PetscSolver0DHandler.
	 *
	 * @param _network The reaction network to use.
	 */
	PetscSolver0DHandler(NetworkType& _network,
		perf::IPerfHandler& _perfHandler, const options::IOptions& options) :
		PetscSolverHandler(_network, _perfHandler, options)
	{
	}

	//! The Destructor
	~PetscSolver0DHandler()
	{
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	createSolverContext(DM& da) override;

	/**
	 * \see ISolverHandler.h
	 */
	void
	initializeSolverContext(DM& da, Mat& J) override;

	/**
	 * \see ISolverHandler.h
	 */
	void
	initializeConcentration(DM& da, Vec& C, DM& oldDA, Vec& oldC) override;

	/**
	 * \see ISolverHandler.h
	 */
	void
	initGBLocation(DM& da, Vec& C) override
	{
		// Doesn't do anything in 0D
		return;
	}

	/**
	 * \see ISolverHandler.h
	 */
	std::vector<
		std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
	getConcVector(DM& da, Vec& C) override;

	/**
	 * \see ISolverHandler.h
	 */
	void
	setConcVector(DM& da, Vec& C,
		std::vector<
			std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>&
			concVector) override;

	/**
	 * \see ISolverHandler.h
	 */
	void
	updateConcentration(TS& ts, Vec& localC, Vec& F, PetscReal ftime) override;

	/**
	 * \see ISolverHandler.h
	 */
	void
	computeJacobian(TS& ts, Vec& localC, Mat& J, PetscReal ftime) override;

	/**
	 * \see ISolverHandler.h
	 */
	IdType
	getSurfacePosition(IdType j = badId, IdType k = badId) const override
	{
		return 0;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setSurfacePosition(IdType pos, IdType j = badId, IdType k = badId) override
	{
		return;
	}

	/**
	.* \see ISolverHandler.h
	 */
	void
	getNetworkTemperature(
		std::vector<double>& temperatures, std::vector<double>& depths) override
	{
		temperatures = temperature;
		depths = std::vector<double>(1, 1.0);
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setSurfaceOffset(int offset) override
	{
		return;
	}
};
// end class PetscSolver0DHandler

} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */
#endif
