#include <xolotl/core/network/LiReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
LiClusterGenerator::LiClusterGenerator(const options::IOptions& options) :
	_maxH(options.getMaxImpurity()),
	_maxV(options.getMaxV())
{
}

LiClusterGenerator::LiClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxH(options.getMaxImpurity()),
	_maxV(options.getMaxV())
{
}
} // namespace network
} // namespace core
} // namespace xolotl
