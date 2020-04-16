#pragma once

#include <Kokkos_Core.hpp>

#include <plsm/Utility.h>

#include <experimental/MemorySpace.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
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

struct ReactionData
{
    using IndexType = ReactionNetworkIndexType;

    ReactionData() = default;

    ReactionData(IndexType numProductionReactions,
            IndexType numDissociationReactions, IndexType numSinkReactions,
            IndexType numReSoReactions, std::size_t numSpeciesNoI,
            IndexType gridSize)
        :
        coeffExtent(numSpeciesNoI + 1),
        numReactions(numProductionReactions + numDissociationReactions + numReSoReactions),
        numRates(numProductionReactions + numDissociationReactions + numSinkReactions + numReSoReactions),
        productionCoeffs("Production Coefficients",
            numProductionReactions, coeffExtent, coeffExtent, 4, coeffExtent),
        dissociationCoeffs("Dissociation Coefficients",
            numDissociationReactions, coeffExtent, 1, 3, coeffExtent),
        resolutionCoeffs("ReSolution Coefficients",
            numReSoReactions, coeffExtent, 1, 3, coeffExtent),
        widths("Reaction Rates", numReactions, numSpeciesNoI),
        rates("Reaction Rates", numRates, gridSize)
    {
    }

    void
    setGridSize(IndexType gridSize) {
        rates = Kokkos::View<double**>("Reaction Rates", numRates, gridSize);
    }

    std::size_t coeffExtent {};
    IndexType numReactions {};
    IndexType numRates {};
    Kokkos::View<double*****> productionCoeffs;
    Kokkos::View<double*****> dissociationCoeffs;
    Kokkos::View<double*****> resolutionCoeffs;
    Kokkos::View<double**> widths;
    Kokkos::View<double**> rates;
    ClusterConnectivity<> connectivity;
};

struct ReactionDataRef
{
    using IndexType = ReactionNetworkIndexType;

    ReactionDataRef() = default;

    KOKKOS_INLINE_FUNCTION
    ReactionDataRef(const ReactionData& data)
        :
        productionCoeffs(data.productionCoeffs),
        dissociationCoeffs(data.dissociationCoeffs),
        resolutionCoeffs(data.resolutionCoeffs),
        widths(data.widths),
        rates(data.rates),
        connectivity(data.connectivity)
    {
    }

    KOKKOS_INLINE_FUNCTION
    auto
    getCoefficients(IndexType reactionId)
    {
        if (reactionId < productionCoeffs.extent(0)) {
            return Kokkos::subview(productionCoeffs, reactionId, Kokkos::ALL,
                Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        }
        else if (reactionId < productionCoeffs.extent(0) + dissociationCoeffs.extent(0)) {
            // TODO: can we use the same indices for dissociation and re-solution as the coefs are the same?
            reactionId -= productionCoeffs.extent(0);
            return Kokkos::subview(dissociationCoeffs, reactionId, Kokkos::ALL,
                Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        }
        else {
            reactionId -= productionCoeffs.extent(0) + dissociationCoeffs.extent(0);
            return Kokkos::subview(resolutionCoeffs, reactionId, Kokkos::ALL,
                Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        }
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

    Kokkos::View<double*****, Kokkos::MemoryUnmanaged> productionCoeffs;
    Kokkos::View<double*****, Kokkos::MemoryUnmanaged> dissociationCoeffs;
    Kokkos::View<double*****, Kokkos::MemoryUnmanaged> resolutionCoeffs;
    Kokkos::View<double**, Kokkos::MemoryUnmanaged> widths;
    Kokkos::View<double**, Kokkos::MemoryUnmanaged> rates;
    ClusterConnectivity<> connectivity;
};
}
}
}
