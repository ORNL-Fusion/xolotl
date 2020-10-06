#pragma once

#include <xolotl/core/temperature/ITemperatureHandler.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
class TemperatureHandler : public ITemperatureHandler
{
public:
	TemperatureHandler() = default;

	virtual ~TemperatureHandler();

	void
	initializeTemperature(int dof,
		network::IReactionNetwork::SparseFillMap& ofillMap,
		network::IReactionNetwork::SparseFillMap& dfillMap) override;

protected:
	/**
	 * The number of degrees of freedom in the network
	 */
	int _dof;
};
} // namespace temperature
} // namespace core
} // namespace xolotl
