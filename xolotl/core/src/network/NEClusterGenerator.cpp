#include <xolotl/core/network/NEReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
NEClusterGenerator::NEClusterGenerator(const options::IOptions& opts) :
	_maxXe(opts.getMaxImpurity()),
	_xeDiffusivity(opts.getXenonDiffusivity()),
	_xeDiffusive(_xeDiffusivity > 0.0),
	_groupingMin(opts.getGroupingMin()),
	_groupingWidth(opts.getGroupingWidthA()),
	_density(opts.getDensity())
{
}

NEClusterGenerator::NEClusterGenerator(
	const options::IOptions& opts, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxXe(opts.getMaxImpurity()),
	_xeDiffusivity(opts.getXenonDiffusivity()),
	_xeDiffusive(_xeDiffusivity > 0.0),
	_groupingMin(opts.getGroupingMin()),
	_groupingWidth(opts.getGroupingWidthA()),
	_density(opts.getDensity())
{
}
} // namespace network
} // namespace core
} // namespace xolotl
