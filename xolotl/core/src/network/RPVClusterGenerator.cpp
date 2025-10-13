#include <xolotl/core/network/RPVReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
RPVClusterGenerator::RPVClusterGenerator(const options::IOptions& options) :
	_maxV(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthV(options.getGroupingWidthA())
{
}

RPVClusterGenerator::RPVClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxV(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthV(options.getGroupingWidthA())
{
}
} // namespace network
} // namespace core
} // namespace xolotl
