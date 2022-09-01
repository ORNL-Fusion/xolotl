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
private:
	//! The position of the surface
	IdType surfacePosition;

public:
	PetscSolver1DHandler() = delete;

	/**
	 * Construct a PetscSolver1DHandler.
	 *
	 * @param _network The reaction network to use.
	 */
	PetscSolver1DHandler(
		NetworkType& _network, const options::IOptions& options) :
		PetscSolverHandler(_network, options),
		surfacePosition(0)
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
	createSolverContext(DM& da);

	/**
	 * \see ISolverHandler.h
	 */
	void
	initializeConcentration(DM& da, Vec& C);

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
		return surfacePosition;
	}

	/**
	 * \see ISolverHandler.h
	 */
	void
	setSurfacePosition(IdType pos, IdType j = -1, IdType k = -1)
	{
		auto oldPos = surfacePosition;
		surfacePosition = pos;
		generateTemperatureGrid(surfacePosition, oldPos);

		return;
	}

	/**
	.* \see ISolverHandler.h
	 */
	void
	getNetworkTemperature(
		std::vector<double>& temperatures, std::vector<double>& depths)
	{
		temperatures = interpolateTemperature(surfacePosition);
		for (auto i = 0; i < temperatures.size(); i++) {
			if (localXS + i == nX + 1)
				depths.push_back(grid[localXS + i] - grid[surfacePosition + 1]);
			else
				depths.push_back(
					(grid[localXS + i + 1] + grid[localXS + i]) / 2.0 -
					grid[surfacePosition + 1]);
		}
	}
};
// end class PetscSolver1DHandler

} /* namespace handler */
} /* namespace solver */
} /* namespace xolotl */
#endif
