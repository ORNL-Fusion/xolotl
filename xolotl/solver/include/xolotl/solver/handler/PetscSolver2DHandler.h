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
	PetscSolver2DHandler(NetworkType& _network,
		perf::IPerfHandler& _perfHandler, const options::IOptions& options) :
		PetscSolverHandler(_network, _perfHandler, options)
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
		return surfacePosition[j];
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setSurfacePosition(IdType pos, IdType j = badId, IdType k = badId) override
	{
		surfacePosition[j] = pos;

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
		for (auto i = 0; i < temperature.size(); i++) {
			if (localXS + i == nX + 1)
				depths.push_back(
					grid[localXS + i] - grid[surfacePosition[localYS] + 1]);
			else
				depths.push_back(
					(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
					grid[surfacePosition[localYS] + 1]);
		}
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
// end class PetscSolver2DHandler

} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */
#endif
