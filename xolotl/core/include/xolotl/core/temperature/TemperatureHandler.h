#ifndef TEMPERATUREHANDLER_H
#define TEMPERATUREHANDLER_H

#include <xolotl/core/temperature/ITemperatureHandler.h>
#include <xolotl/factory/temperature/TemperatureHandlerFactory.h>
#include <xolotl/options/Options.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
class TemperatureHandler : public ITemperatureHandler
{
public:
	TemperatureHandler(const options::Options& options);

	virtual ~TemperatureHandler()
	{
	}

	virtual void
	initializeTemperature(const int dof,
		network::IReactionNetwork::SparseFillMap& ofillMap,
		network::IReactionNetwork::SparseFillMap& dfillMap) override
	{
		_strategy->initializeTemperature(dof, ofillMap, dfillMap);
	}

	virtual double
	getTemperature(const plsm::SpaceVector<double, 3>& fraction,
		double currentTime) const override
	{
		return _strategy->getTemperature(fraction, currentTime);
	}

	virtual void
	setTemperature(double* solution) override
	{
		_strategy->setTemperature(solution);
	}

	virtual void
	setHeatCoefficient(double coef) override
	{
		_strategy->setHeatCoefficient(coef);
	}

	virtual void
	setHeatConductivity(double cond) override
	{
		_strategy->setHeatConductivity(cond);
	}

	virtual void
	updateSurfacePosition(int surfacePos) override
	{
		_strategy->updateSurfacePosition(surfacePos);
	}

	virtual void
	computeTemperature(double** concVector, double* updatedConcOffset,
		double hxLeft, double hxRight, int xi, double sy = 0.0, int iy = 0,
		double sz = 0.0, int iz = 0) override
	{
		_strategy->computeTemperature(
			concVector, updatedConcOffset, hxLeft, hxRight, xi, sy, iy, sz, iz);
	}

	virtual bool
	computePartialsForTemperature(double* val, int* indices, double hxLeft,
		double hxRight, int xi, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) override
	{
		return _strategy->computePartialsForTemperature(
			val, indices, hxLeft, hxRight, xi, sy, iy, sz, iz);
	}

private:
	std::shared_ptr<ITemperatureHandler> _strategy;
};
} // namespace temperature
} // namespace core
} // namespace xolotl

#endif
