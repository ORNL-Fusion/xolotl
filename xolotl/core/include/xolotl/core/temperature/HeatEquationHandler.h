#pragma once

#include <xolotl/core/temperature/TemperatureHandler.h>
#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
class HeatEquationHandler : public TemperatureHandler
{
public:
	HeatEquationHandler() = delete;

	HeatEquationHandler(double flux, double bulkTemp, int dim);

	HeatEquationHandler(const options::IOptions& options);

	virtual ~HeatEquationHandler();

	double
	getTemperature(
		const plsm::SpaceVector<double, 3>&, double time) const override;

	void
	setTemperature(double* solution) override
	{
		localTemperature = solution[this->_dof];
	}

	void
	setHeatCoefficient(double coef) override
	{
		heatCoef = coef;
	}

	void
	setHeatConductivity(double cond) override
	{
		heatConductivity = cond;
	}

	void
	updateSurfacePosition(int surfacePos) override
	{
		surfacePosition = surfacePos;
	}

	void
	computeTemperature(double** concVector, double* updatedConcOffset,
		double hxLeft, double hxRight, int xi, double sy = 0.0, int iy = 0,
		double sz = 0.0, int iz = 0) override;

	bool
	computePartialsForTemperature(double* val, int* indices, double hxLeft,
		double hxRight, int xi, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) override;

private:
	/**
	 * The heat flux in W.m-2
	 */
	double heatFlux;

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
	 * The heat coefficient
	 */
	double heatCoef;

	/**
	 * The heat conductivity
	 */
	double heatConductivity;

	/**
	 * Zero flux enables constant temperature behavior
	 */
	bool zeroFlux{false};

	/**
	 * Number of dimensions for the current simulation
	 */
	int dimension;

	/**
	 * Hang on to single allocation for use in computeTemperature()
	 */
	std::vector<std::array<double, 2>> oldConcBox;
};
} // namespace temperature
} // namespace core
} // namespace xolotl
