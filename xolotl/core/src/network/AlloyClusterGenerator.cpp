#include <xolotl/core/network/AlloyReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
AlloyClusterGenerator::AlloyClusterGenerator(const options::IOptions& options) :
	_maxV(options.getMaxV()),
	_maxI(options.getMaxI()),
	_maxSize(options.getMaxImpurity()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidth(options.getGroupingWidthA())
{
}

AlloyClusterGenerator::AlloyClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxV(options.getMaxV()),
	_maxI(options.getMaxI()),
	_maxSize(options.getMaxImpurity()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidth(options.getGroupingWidthA())
{
}
} // namespace network
} // namespace core
} // namespace xolotl
