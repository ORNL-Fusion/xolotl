#ifndef PETSCSOLVER1DHANDLER_H
#define PETSCSOLVER1DHANDLER_H

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
 * needed to solve the DR equations in 1D using PETSc from Argonne National
 * Laboratory.
 */
class PetscSolver1DHandler : public PetscSolverHandler
{
public:
	PetscSolver1DHandler() = delete;

	/**
	 * Construct a PetscSolver1DHandler.
	 *
	 * @param _network The reaction network to use.
	 */
	PetscSolver1DHandler(NetworkType& _network,
		perf::IPerfHandler& _perfHandler, const options::IOptions& options) :
		PetscSolverHandler(_network, _perfHandler, options)
	{
	}

	//! The Destructor
	~PetscSolver1DHandler()
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
	initGBLocation(DM& da, Vec& C) override;

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
		temperatures = interpolateTemperature();
		for (auto i = 0; i < temperatures.size(); i++) {
			if (localXS + i == nX + 1)
				depths.push_back(grid[localXS + i] - grid[1]);
			else
				depths.push_back(
					(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
					grid[1]);
		}
	}
};
// end class PetscSolver1DHandler

} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */
#endif
