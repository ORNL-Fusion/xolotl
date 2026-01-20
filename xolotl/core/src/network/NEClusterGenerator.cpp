#include <xolotl/core/network/NEReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
NEClusterGenerator::NEClusterGenerator(const options::IOptions& options) :
	_maxXe(options.getMaxImpurity()),
	_maxV(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthXe(options.getGroupingWidthA()),
	_groupingWidthV(options.getGroupingWidthB()),
	_xeDiffusivity(options.getXenonDiffusivity()),
	_xeDiffusive(_xeDiffusivity > 0.0),
	_density(options.getDensity())
{
}

NEClusterGenerator::NEClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxXe(options.getMaxImpurity()),
	_maxV(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthXe(options.getGroupingWidthA()),
	_groupingWidthV(options.getGroupingWidthB()),
	_xeDiffusivity(options.getXenonDiffusivity()),
	_xeDiffusive(_xeDiffusivity > 0.0),
	_density(options.getDensity())
{
}
} // namespace network
} // namespace core
} // namespace xolotl
