#include <xolotl/core/network/FeCrReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
FeCrClusterGenerator::FeCrClusterGenerator(const options::IOptions& options) :
	_minJunction(options.getMaxV()),
	_maxI(options.getMaxI()),
	_maxSize(options.getMaxD()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidth(options.getGroupingWidthA()),
	_maxHe(1)
{
}

FeCrClusterGenerator::FeCrClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_minJunction(options.getMaxV()),
	_maxI(options.getMaxI()),
	_maxSize(options.getMaxD()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidth(options.getGroupingWidthA()),
	_maxHe(1)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
