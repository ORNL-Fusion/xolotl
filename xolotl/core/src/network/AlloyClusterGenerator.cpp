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
	_maxLoopSize(options.getMaxImpurity()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthA(options.getGroupingWidthA()),
	_groupingWidthB(options.getGroupingWidthB()),
	_hevRatio(options.getHeVRatio())
{
	// Check if SSBM is used
	auto process = options.getProcesses();
	if (process["largeBubble"]) {
		_maxLoopSize = 2000;
	}
}

AlloyClusterGenerator::AlloyClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxV(options.getMaxV()),
	_maxI(options.getMaxI()),
	_maxSize(options.getMaxImpurity()),
	_maxLoopSize(options.getMaxImpurity()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthA(options.getGroupingWidthA()),
	_groupingWidthB(options.getGroupingWidthB()),
	_hevRatio(options.getHeVRatio())
{
	// Check if SSBM is used
	auto process = options.getProcesses();
	if (process["largeBubble"]) {
		_maxLoopSize = 2000;
	}
}
} // namespace network
} // namespace core
} // namespace xolotl
