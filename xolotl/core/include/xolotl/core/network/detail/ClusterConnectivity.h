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

	ClusterConnectivity&
	operator=(const ClusterConnectivity& other)
	{
		Crs::operator=(other);
		if (other._tr.row_map.size() == 0 && this->row_map.size() > 0) {
			_tr = Crs();
			Kokkos::transpose_crs(_tr, *this);
			_trEntries = typename Crs::entries_type(
				"transposed entries map", this->entries.size());
			auto conn = *this;
			auto nRows = conn.row_map.size() - 1;
			Kokkos::parallel_for(
				"ClusterConnectivity::assignTransposeEntries", nRows,
				KOKKOS_LAMBDA(IndexType i) {
					auto pBegin = conn.row_map(i);
					auto pEnd = conn.row_map(i + 1);
					for (auto p = pBegin; p < pEnd; ++p) {
						auto j = conn.entries(p);
						conn._trEntries(conn.getPosition(j, i, conn._tr)) = p;
					}
				});
			Kokkos::fence();
			_avgRowSize = static_cast<double>(conn.entries.size()) /
					static_cast<double>(nRows) +
				0.5;
		}
		else {
			_tr = other._tr;
			_trEntries = other._trEntries;
			_avgRowSize = other._avgRowSize;
		}
		return *this;
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
		if (getRowSize(rowId) > _avgRowSize) {
			auto trPos = getPosition(columnId, rowId, _tr);
			if (trPos == invalidNetworkIndex) {
				return invalidNetworkIndex;
			}
			return _trEntries(trPos);
		}
		else {
			return getPosition(rowId, columnId, *this);
		}
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		std::uint64_t ret{};

		ret += this->row_map.required_allocation_size(this->row_map.extent(0));
		ret += this->entries.required_allocation_size(this->entries.extent(0));
		ret += _tr.row_map.required_allocation_size(_tr.row_map.extent(0));
		ret += _tr.entries.required_allocation_size(_tr.entries.extent(0));
		ret += _trEntries.required_allocation_size(_trEntries.extent(0));
		ret += sizeof(_avgRowSize);

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

private:
	Crs _tr;
	typename Crs::entries_type _trEntries;
	IndexType _avgRowSize{};
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
