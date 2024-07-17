#include <xolotl/core/network/T91ReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
T91ClusterGenerator::T91ClusterGenerator(const options::IOptions& options) :
	_maxHe(options.getMaxImpurity()),
	_maxV(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthHe(options.getGroupingWidthA()),
	_groupingWidthV(options.getGroupingWidthB())
{
}

T91ClusterGenerator::T91ClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxHe(options.getMaxImpurity()),
	_maxV(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthHe(options.getGroupingWidthA()),
	_groupingWidthV(options.getGroupingWidthB())
{
}
} // namespace network
} // namespace core
} // namespace xolotl
