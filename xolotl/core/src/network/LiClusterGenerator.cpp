#include <xolotl/core/network/LiReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
LiClusterGenerator::LiClusterGenerator(const options::IOptions& options) :
	_maxH(options.getMaxImpurity())
{
}

LiClusterGenerator::LiClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxH(options.getMaxImpurity())
{
}
} // namespace network
} // namespace core
} // namespace xolotl
