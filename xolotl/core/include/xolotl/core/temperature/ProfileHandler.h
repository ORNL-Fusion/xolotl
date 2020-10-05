#pragma once

#include <string>

#include <xolotl/core/temperature/TemperatureHandler.h>
#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of a temperature changing with time.
 */
class ProfileHandler : public TemperatureHandler
{
private:
	/**
	 * The name of the file were the profile is stored.
	 */
	std::string tempFile;

	/**
	 * Vector to hold the time read from the input
	 * temperature file.
	 */
	std::vector<double> time;

	/**
	 * Vector to hold the temperature read from the input
	 * temperature file.
	 */
	std::vector<double> temp;

public:
	ProfileHandler() = delete;

	/**
	 * The constructor.
	 *
	 * @param profileFileName The name of the profile file
	 */
	ProfileHandler(const std::string& profileFileName);

    /**
     * Construct from options
     */
    ProfileHandler(const options::IOptions& options);

	/**
	 * The destructor.
	 */
	virtual ~ProfileHandler();

	/**
	 * This operation initializes the ofill and dfill arrays so that the
	 * temperature is connected correctly in the solver.
	 * It also reads in the time and temperature data from the input
	 * temperature file that was specified by the command line.
	 *
	 * \see ITemperatureHandler.h
	 */
	void
	initializeTemperature(const int dof,
		network::IReactionNetwork::SparseFillMap& ofillMap,
		network::IReactionNetwork::SparseFillMap& dfillMap) override;

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 * It linearly interpolates the data read from the input
	 * temperature file.
	 *
	 * \see ITemperatureHandler.h
	 */
	double
	getTemperature(
		const plsm::SpaceVector<double, 3>&, double currentTime) const override;

	/**
	 * This operation sets the temperature given by the solver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	void
	setTemperature(double* solution) override
	{
		return;
	}

	/**
	 * This operation sets the heat coefficient to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	void
	setHeatCoefficient(double coef) override
	{
		return;
	}

	/**
	 * This operation sets the heat conductivity to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	void
	setHeatConductivity(double cond) override
	{
		return;
	}

	/**
	 * This operation sets the surface position.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	void
	updateSurfacePosition(int surfacePos) override
	{
		return;
	}

	/**
	 * Compute the flux due to the heat equation.
	 * This method is called by the RHSFunction from the PetscSolver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	void
	computeTemperature(double** concVector, double* updatedConcOffset,
		double hxLeft, double hxRight, int xi, double sy = 0.0, int iy = 0,
		double sz = 0.0, int iz = 0) override
	{
		return;
	}

	/**
	 * Compute the partials due to the heat equation.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	bool
	computePartialsForTemperature(double* val, int* indices, double hxLeft,
		double hxRight, int xi, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) override
	{
		return false;
	}
};
// end class ProfileHandler

} // namespace temperature
} // namespace core
} // namespace xolotl
