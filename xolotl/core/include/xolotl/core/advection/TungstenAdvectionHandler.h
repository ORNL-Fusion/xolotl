#pragma once

#include <array>

#include <xolotl/core/advection/SurfaceAdvectionHandler.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace advection
{
/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 */
class TungstenAdvectionHandler : public SurfaceAdvectionHandler
{
public:
	TungstenAdvectionHandler(const std::array<double, 7>& sinkStrength);

	virtual ~TungstenAdvectionHandler();

	/**
	 * This function initializes the list of clusters that will move through
	 * advection for a (100) tungsten material.
	 *
	 * @param network The network
	 * @param ofillMap Map of connectivity for advecting clusters.
	 * of the advecting clusters
	 */
	void
	initialize(network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& ofillMap) override;

protected:
	std::array<double, 7> _sinkStrength{};
};
} // namespace advection
} // namespace core
} // namespace xolotl
