#pragma once

#include <Kokkos_Core.hpp>

#include <plsm/Utility.h>

#include <xolotl/core/reactants/MemorySpace.h>
#include <xolotl/core/reactants/ReactionNetworkTraits.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
using CoefficientsView = Kokkos::View<double*****>;
using CoefficientsViewUnmanaged =
    Kokkos::View<double*****, Kokkos::MemoryUnmanaged>;

template <typename TMemSpace = DefaultMemorySpace>
class ClusterConnectivity :
    public Kokkos::Crs<ReactionNetworkIndexType, TMemSpace>
{
public:
    using IndexType = ReactionNetworkIndexType;

private:
    static constexpr IndexType invalidIndex = InvalidIndex::value;

public:
    using Crs = Kokkos::Crs<IndexType, TMemSpace>;

    using HostMirror =
        ClusterConnectivity<
            typename Kokkos::View<int*, TMemSpace>::traits::host_mirror_space>;

    ClusterConnectivity&
    operator=(const ClusterConnectivity& other)
    {
        Crs::operator=(other);
        if (other._tr.row_map.size() == 0 && this->row_map.size() > 0) {
            _tr = Crs();
            Kokkos::transpose_crs(_tr, *this);
            _trEntries = typename Crs::entries_type("transposed entries map",
                this->entries.size());
            auto conn = *this;
            auto nRows = conn.row_map.size() - 1;
            Kokkos::parallel_for(nRows, KOKKOS_LAMBDA (IndexType i) {
                auto pBegin = conn.row_map(i);
                auto pEnd = conn.row_map(i+1);
                for (auto p = pBegin; p < pEnd; ++p) {
                    auto j = conn.entries(p);
                    conn._trEntries(conn.getPosition(j, i, conn._tr)) = p;
                }
            });
            Kokkos::fence();
            _avgRowSize = static_cast<double>(conn.entries.size()) /
                static_cast<double>(nRows) + 0.5;
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
        if (rowId == invalidIndex || columnId == invalidIndex) {
            return;
        }
        if (this->entries.size() == 0) {
            Kokkos::atomic_increment(&this->row_map(rowId));
        }
        else {
            for (auto id = this->row_map(rowId);
                    !Kokkos::atomic_compare_exchange_strong(&this->entries(id),
                        invalidIndex, columnId); ++id) {
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
            if (trPos == invalidIndex) {
                return invalidIndex;
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
        std::uint64_t ret {};

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
        for (auto pos = crs.row_map(rowId); pos < crs.row_map(rowId+1); ++pos) {
            if (crs.entries(pos) == columnId) {
                return pos;
            }
        }
        return invalidIndex;
    }

    KOKKOS_INLINE_FUNCTION
    IndexType
    getRowSize(IndexType rowId) const
    {
        return this->row_map(rowId+1) - this->row_map(rowId);
    }

private:
    Crs _tr;
    typename Crs::entries_type _trEntries;
    IndexType _avgRowSize {};
};


template <typename TNetwork>
struct ReactionData
{
    using NetworkType = TNetwork;
    using IndexType = ReactionNetworkIndexType;
    using ReactionTypes = ReactionTypeList<NetworkType>;

    static constexpr std::size_t numReactionTypes =
        std::tuple_size<ReactionTypes>::value;

    static constexpr std::size_t numSpeciesNoI =
        NetworkType::getNumberOfSpeciesNoI();

    ReactionData() = default;

    ReactionData(IndexType nReactions, IndexType gridSize,
            const Kokkos::Array<IndexType, numReactionTypes+1>& rBeginIds)
        :
        numReactions(nReactions),
        widths("Reaction Widths", numReactions, numSpeciesNoI),
        rates("Reaction Rates", numReactions, gridSize),
        reactionBeginIndices(rBeginIds)
    {
    }

    std::uint64_t
    getDeviceMemorySize() const noexcept
    {
        std::uint64_t ret = sizeof(numReactions);
        ret += sizeof(reactionBeginIndices);
        ret += widths.required_allocation_size(widths.extent(0),
            widths.extent(1));
        ret += rates.required_allocation_size(rates.extent(0), rates.extent(1));

        for (std::size_t r = 0; r < numReactionTypes; ++r) {
            ret += coeffs[r].required_allocation_size(coeffs[r].extent(0),
                coeffs[r].extent(1), coeffs[r].extent(2),
                coeffs[r].extent(3), coeffs[r].extent(4));
        }

        ret += connectivity.getDeviceMemorySize();

        return ret;
    }

    void
    setGridSize(IndexType gridSize)
    {
        rates = Kokkos::View<double**>("Reaction Rates", numReactions, gridSize);
    }

    IndexType numReactions {};
    Kokkos::View<double**> widths;
    Kokkos::View<double**> rates;
    Kokkos::Array<IndexType, numReactionTypes+1> reactionBeginIndices;
    Kokkos::Array<CoefficientsView, numReactionTypes> coeffs;
    ClusterConnectivity<> connectivity;
};

template <typename TNetwork>
struct ReactionDataRef
{
    using NetworkType = TNetwork;
    using IndexType = ReactionNetworkIndexType;
    using ReactionTypes = ReactionTypeList<NetworkType>;

    static constexpr std::size_t numReactionTypes =
        std::tuple_size<ReactionTypes>::value;

    ReactionDataRef() = default;

    KOKKOS_INLINE_FUNCTION
    ReactionDataRef(const ReactionData<NetworkType>& data)
        :
        widths(data.widths),
        rates(data.rates),
        reactionBeginIndices(data.reactionBeginIndices),
        connectivity(data.connectivity)
    {
        for (std::size_t r = 0; r < numReactionTypes; ++r) {
            coeffs[r] = data.coeffs[r];
        }
    }

    KOKKOS_INLINE_FUNCTION
    auto
    getCoefficients(IndexType reactionId)
    {
        std::size_t r = 0;
        for (; r < numReactionTypes; ++r) {
            if (reactionId < reactionBeginIndices[r+1]) {
                break;
            }
        }
        assert(r < numReactionTypes);
        return Kokkos::subview(coeffs[r], reactionId - reactionBeginIndices[r],
            Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    }

    KOKKOS_INLINE_FUNCTION
    auto
    getWidths(IndexType reactionId)
    {
        return Kokkos::subview(widths, reactionId, Kokkos::ALL);
    }

    KOKKOS_INLINE_FUNCTION
    auto
    getRates(IndexType reactionId)
    {
        return Kokkos::subview(rates, reactionId, Kokkos::ALL);
    }

    Kokkos::View<double**, Kokkos::MemoryUnmanaged> widths;
    Kokkos::View<double**, Kokkos::MemoryUnmanaged> rates;
    Kokkos::Array<IndexType, numReactionTypes+1> reactionBeginIndices;
    Kokkos::Array<CoefficientsViewUnmanaged, numReactionTypes> coeffs;
    //TODO: Enable unmanaged version of connectivity
    ClusterConnectivity<> connectivity;
};
}
}
}
