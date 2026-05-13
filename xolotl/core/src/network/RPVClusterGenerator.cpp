#include <xolotl/core/network/RPVReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
RPVClusterGenerator::RPVClusterGenerator(const options::IOptions& options) :
	_maxI(options.getMaxI()),
	_maxSize(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidth(options.getGroupingWidthA())
{
}

RPVClusterGenerator::RPVClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxI(options.getMaxI()),
	_maxSize(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidth(options.getGroupingWidthA())
{
}
} // namespace network
} // namespace core
} // namespace xolotl
