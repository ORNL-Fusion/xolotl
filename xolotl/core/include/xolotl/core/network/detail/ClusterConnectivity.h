#pragma once

#include <Kokkos_Core.hpp>

#include <xolotl/core/network/ReactionNetworkTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
/**
 * @brief Connectivity object that stores which
 * cluster interacts with which, needed to define the sparse
 * Jacobian. Relies on Kokkos Compressed Row Storage array.
 *
 * @tparam TMemSpace The memory layout type
 */
template <typename TMemSpace = plsm::DefaultMemSpace>
class ClusterConnectivity :
	public Kokkos::Crs<ReactionNetworkIndexType, TMemSpace>
{
public:
	using IndexType = ReactionNetworkIndexType;

public:
	using Crs = Kokkos::Crs<IndexType, TMemSpace>;

	using HostMirror = ClusterConnectivity<
		typename Kokkos::View<int*, TMemSpace>::traits::host_mirror_space>;

	KOKKOS_INLINE_FUNCTION
	IndexType
	getNumberOfRows() const
	{
		return this->row_map.size() - 1;
	}

	KOKKOS_INLINE_FUNCTION
	void
	add(IndexType rowId, IndexType columnId) const
	{
		if (rowId == invalidNetworkIndex || columnId == invalidNetworkIndex) {
			return;
		}
		if (this->entries.size() == 0) {
			Kokkos::atomic_increment(&this->row_map(rowId));
		}
		else {
			for (auto id = this->row_map(rowId);
				 !Kokkos::atomic_compare_exchange_strong(
					 &this->entries(id), invalidNetworkIndex, columnId);
				 ++id) {
				if (this->entries(id) == columnId) {
					break;
				}
			}
		}
	}

	KOKKOS_INLINE_FUNCTION
	IndexType
	operator()(IndexType rowId, IndexType columnId) const
	{
        return getPosition(rowId, columnId, *this);
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		std::uint64_t ret{};

		ret += this->row_map.required_allocation_size(this->row_map.extent(0));
		ret += this->entries.required_allocation_size(this->entries.extent(0));

		return ret;
	}

private:
	KOKKOS_INLINE_FUNCTION
	IndexType
	getPosition(IndexType rowId, IndexType columnId, const Crs& crs) const
	{
		for (auto pos = crs.row_map(rowId); pos < crs.row_map(rowId + 1);
			 ++pos) {
			if (crs.entries(pos) == columnId) {
				return pos;
			}
		}
		return invalidNetworkIndex;
	}

	KOKKOS_INLINE_FUNCTION
	IndexType
	getRowSize(IndexType rowId) const
	{
		return this->row_map(rowId + 1) - this->row_map(rowId);
	}
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
