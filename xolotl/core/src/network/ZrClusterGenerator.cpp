#include <xolotl/core/network/ZrReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
ZrClusterGenerator::ZrClusterGenerator(const options::IOptions& options) :
	_maxV(options.getMaxV()),
	_maxB(options.getMaxImpurity()),
	_maxI(options.getMaxI()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidth(options.getGroupingWidthA()),
	_transitionSize(options.getTransitionSize())
{
}

ZrClusterGenerator::ZrClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxV(options.getMaxV()),
	_maxB(options.getMaxImpurity()),
	_maxI(options.getMaxI()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidth(options.getGroupingWidthA()),
	_transitionSize(options.getTransitionSize())
{
}
} // namespace network
} // namespace core
} // namespace xolotl
