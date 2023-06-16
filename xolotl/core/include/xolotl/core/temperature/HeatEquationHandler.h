#pragma once

#include <xolotl/core/temperature/TemperatureHandler.h>
#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
/**
 * This class realizes the TemperatureHandler, it is responsible for the
 * handling of a temperature following the heat equation.
 */
class HeatEquationHandler : public TemperatureHandler
{
public:
	HeatEquationHandler() = delete;

	/**
	 * The constructor.
	 *
	 * @param flux The heat flux
	 * @param bulkTemp The bulk temperature
	 * @param dim the number of dimensions for the simulation
	 * @param filename The name of the file containing the heat flux time
	 * profile
	 */
	HeatEquationHandler(
		double flux, double bulkTemp, int dim, std::string filename = "");

	/**
	 * Construct from options
	 */
	HeatEquationHandler(const options::IOptions& options);

	/**
	 * The destructor.
	 */
	virtual ~HeatEquationHandler();

	/**
	 * \see ITemperatureHandler.h
	 */
	double
	getTemperature(
		const plsm::SpaceVector<double, 3>&, double time) const override;

	/**
	 * \see ITemperatureHandler.h
	 */
	void
	setTemperature(double* solution) override
	{
		localTemperature = solution[this->_dof];
	}

	/**
	 * \see ITemperatureHandler.h
	 */
	void
	setHeatCoefficient(double coef) override
	{
		heatCoef = coef;
	}

	/**
	 * \see ITemperatureHandler.h
	 */
	void
	setHeatConductivity(double cond) override
	{
		heatConductivity = cond;
	}

	/**
	 * \see ITemperatureHandler.h
	 */
	void
	updateSurfacePosition(int surfacePos, std::vector<double> grid) override
	{
		surfacePosition = surfacePos;
		// keep the grid
		xGrid = grid;
		bulkPosition = grid.size() - 3;
	}

	/**
	 * \see ITemperatureHandler.h
	 */
	void
	computeTemperature(double currentTime, double** concVector,
		double* updatedConcOffset, double hxLeft, double hxRight, int xi,
		double sy = 0.0, int iy = 0, double sz = 0.0, int iz = 0) override;

	/**
	 * \see ITemperatureHandler.h
	 */
	bool
	computePartialsForTemperature(double currentTime, double** concVector,
		double* val, IdType* indices, double hxLeft, double hxRight, int xi,
		double sy = 0.0, int iy = 0, double sz = 0.0, int iz = 0) override;

	/**
	 * Get the heat flux at this time.
	 *
	 * @param currentTime The current time
	 * @return The heat flux
	 */
	double
	getHeatFlux(double currentTime) override;

private:
	/**
	 * The bulk temperature in Kelvin
	 */
	double bulkTemperature;

	/**
	 * The local temperature in Kelvin
	 */
	double localTemperature;

	/**
	 * The surface position
	 */
	int surfacePosition;

	/**
	 * The bulk position
	 */
	int bulkPosition;

	/**
	 * The heat coefficient
	 */
	double heatCoef;

	/**
	 * The heat conductivity
	 */
	double heatConductivity;

	/**
	 * The interface location
	 */
	double interfaceLoc;

	/**
	 * Zero flux enables constant temperature behavior
	 */
	bool zeroFlux{false};

	/**
	 * ELM flux enables time-dependent heat flux
	 */
	bool elmFlux{false};

	/**
	 * Number of dimensions for the current simulation
	 */
	int dimension;
	/**
	 * The name of the file were the profile is stored.
	 */
	std::string fluxFile;

	/**
	 * Vector to hold the time read from the input
	 * temperature file.
	 */
	std::vector<double> time;

	/**
	 * Vector to hold the temperature read from the input
	 * temperature file.
	 */
	std::vector<double> flux;

	/**
	 * Heat conductivity fit
	 */
	double A = 10.846, B = -184.22, C = 872.47;

	/**
	 * Hang on to single allocation for use in computeTemperature()
	 */
	std::vector<std::array<double, 2>> oldConcBox;

	/**
	 * Get the spatially dependent part of the heat conductivity.
	 *
	 * @param xi The grid index
	 * @return Alpha
	 */
	double
	getLocalHeatAlpha(int xi) const;

	/**
	 * Get the temperature dependent part of the heat conductivity.
	 *
	 * @param temp The temperature
	 * @return Beta
	 */
	double
	getLocalHeatBeta(double temp) const;

	/**
	 * Get the inverse of the temperature dependent heat capacity times density
	 * (1.0/(\rho C_v)).
	 *
	 * @param temp The temperature
	 * @return Gamma
	 */
	double
	getLocalHeatGamma(double temp) const;

	/**
	 * Get the first temperature derivative of Beta.
	 *
	 * @param temp The temperature
	 * @return The derivative
	 */
	double
	getDBeta(double temp) const;

	/**
	 * Get the second temperature derivative of Beta.
	 *
	 * @param temp The temperature
	 * @return The second derivative
	 */
	double
	getDDBeta(double temp) const;

	/**
	 * Get the spatial derivative of Alpha.
	 *
	 * @param xi The grid index
	 * @return The derivative
	 */
	double
	getDAlpha(int xi) const;

	/**
	 * Get the temperature derivative of Gamma.
	 *
	 * @param temp The temperature
	 * @return The derivative
	 */
	double
	getDGamma(double temp) const;

	/**
	 * Get the bulk heat flux.
	 *
	 * @param temp The temperature
	 * @return The flux
	 */
	double
	getBulkHeatFlux(double temp) const;

	/**
	 * Get the bulk heat flux derivative.
	 *
	 * @param temp The temperature
	 */
	double
	getBulkHeatFluxDerivative(double temp) const;
};
} // namespace temperature
} // namespace core
} // namespace xolotl
