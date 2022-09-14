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
	 * @param dim the number of dimensions fo the simulation
	 */
	HeatEquationHandler(double flux, double bulkTemp, int dim);

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

	void
	setTemperature(Kokkos::View<const double*> solution) override;

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
	updateSurfacePosition(int surfacePos) override
	{
		surfacePosition = surfacePos;
	}

	/**
	 * \see ITemperatureHandler.h
	 */
	void
	computeTemperature(double** concVector, double* updatedConcOffset,
		double hxLeft, double hxRight, int xi, double sy = 0.0, int iy = 0,
		double sz = 0.0, int iz = 0) override;

	void
	computeTemperature(Kokkos::View<const double*>* concVector,
		Kokkos::View<double*> updatedConcOffset, double hxLeft, double hxRight,
		int xi, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) override;

	/**
	 * \see ITemperatureHandler.h
	 */
	bool
	computePartialsForTemperature(double* val, IdType* indices, double hxLeft,
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
