#pragma once

#include <xolotl/core/network/IReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
/**
 * @brief Structure used by the reactions to store
 * reactant and product ids.
 */
struct ClusterSet
{
	using IndexType = IReactionNetwork::IndexType;
	static constexpr IndexType invalidIndex = IReactionNetwork::invalidIndex();

	IndexType cluster0{invalidIndex};
	IndexType cluster1{invalidIndex};
	IndexType cluster2{invalidIndex};
	IndexType cluster3{invalidIndex};

	ClusterSet() = default;

	KOKKOS_INLINE_FUNCTION
	ClusterSet(IndexType cl0, IndexType cl1, IndexType cl2 = invalidIndex,
		IndexType cl3 = invalidIndex) :
		cluster0{cl0},
		cluster1{cl1},
		cluster2{cl2},
		cluster3{cl3}
	{
	}

	KOKKOS_INLINE_FUNCTION
	bool
	valid() const noexcept
	{
		return cluster0 != invalidIndex;
	}
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
