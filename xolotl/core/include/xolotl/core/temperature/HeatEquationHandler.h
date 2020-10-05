#pragma once

#include <xolotl/core/temperature/ITemperatureHandler.h>
#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
class HeatEquationHandler : public ITemperatureHandler
{
public:
	HeatEquationHandler() = delete;

	HeatEquationHandler(double flux, double bulkTemp, int dim);

	HeatEquationHandler(const options::IOptions& options);

	virtual ~HeatEquationHandler();

	void
	initializeTemperature(const int _dof,
		network::IReactionNetwork::SparseFillMap& ofillMap,
		network::IReactionNetwork::SparseFillMap& dfillMap) override
	{
		// Set dof
		dof = _dof;

		// Add the temperature to ofill
		ofillMap[dof].emplace_back(dof);

		// Add the temperature to dfill
		dfillMap[dof].emplace_back(dof);

		return;
	}

	double
	getTemperature(
		const plsm::SpaceVector<double, 3>&, double time) const override;

	void
	setTemperature(double* solution) override
	{
		localTemperature = solution[dof];
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
	 * The number of degrees of freedom in the network
	 */
	int dof;

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

	bool zeroFlux{false};
	int dimension;
	std::vector<std::array<double, 2>> oldConcBox;
};
} // namespace temperature
} // namespace core
} // namespace xolotl
