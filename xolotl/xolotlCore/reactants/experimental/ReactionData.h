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
class ClusterConnectivity : public Kokkos::Crs<std::size_t, TMemSpace>
{
    static constexpr std::size_t invalid = plsm::invalid<std::size_t>;
public:
    using Crs = Kokkos::Crs<std::size_t, TMemSpace>;

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
            Kokkos::parallel_for(nRows, KOKKOS_LAMBDA (std::size_t i) {
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
    add(std::size_t rowId, std::size_t columnId) const
    {
        if (rowId == invalid || columnId == invalid) {
            return;
        }
        if (this->entries.size() == 0) {
            Kokkos::atomic_increment(&this->row_map(rowId));
        }
        else {
            for (auto id = this->row_map(rowId);
                    !Kokkos::atomic_compare_exchange_strong(&this->entries(id),
                        invalid, columnId); ++id) {
                if (this->entries(id) == columnId) {
                    break;
                }
            }
        }
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t
    operator()(std::size_t rowId, std::size_t columnId) const
    {
        if (getRowSize(rowId) > _avgRowSize) {
            auto trPos = getPosition(columnId, rowId, _tr);
            if (trPos == invalid) {
                return invalid;
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
    std::size_t
    getPosition(std::size_t rowId, std::size_t columnId, const Crs& crs) const
    {
        for (auto pos = crs.row_map(rowId); pos < crs.row_map(rowId+1); ++pos) {
            if (crs.entries(pos) == columnId) {
                return pos;
            }
        }
        return invalid;
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t
    getRowSize(std::size_t rowId) const
    {
        return this->row_map(rowId+1) - this->row_map(rowId);
    }

private:
    Crs _tr;
    typename Crs::entries_type _trEntries;
    std::size_t _avgRowSize {};
};

struct ReactionData
{
    ReactionData() = default;

    ReactionData(std::size_t numProductionReactions,
            std::size_t numDissociationReactions, std::size_t numSpeciesNoI,
            std::size_t gridSize)
        :
        coeffExtent(numSpeciesNoI + 1),
        numReactions(numProductionReactions + numDissociationReactions),
        productionCoeffs("Production Coefficients",
            numProductionReactions, coeffExtent, coeffExtent, 4, coeffExtent),
        dissociationCoeffs("Dissociation Coefficients",
            numDissociationReactions, coeffExtent, 1, 3, coeffExtent),
        rates("Reaction Rates", numReactions, gridSize)
    {
    }

    void
    setGridSize(std::size_t gridSize) {
        rates = Kokkos::View<double**>("Reaction Rates", numReactions, gridSize);
    }

    std::size_t coeffExtent {};
    std::size_t numReactions {};
    Kokkos::View<double*****> productionCoeffs;
    Kokkos::View<double*****> dissociationCoeffs;
    Kokkos::View<double**> rates;
    ClusterConnectivity<> connectivity;
};

struct ReactionDataRef
{
    ReactionDataRef() = default;

    KOKKOS_INLINE_FUNCTION
    ReactionDataRef(const ReactionData& data)
        :
        productionCoeffs(data.productionCoeffs),
        dissociationCoeffs(data.dissociationCoeffs),
        rates(data.rates),
        connectivity(data.connectivity)
    {
    }

    KOKKOS_INLINE_FUNCTION
    auto
    getCoefficients(std::size_t reactionId)
    {
        if (reactionId < productionCoeffs.extent(0)) {
            return Kokkos::subview(productionCoeffs, reactionId, Kokkos::ALL,
                Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        }
        else {
            reactionId -= productionCoeffs.extent(0);
            return Kokkos::subview(dissociationCoeffs, reactionId, Kokkos::ALL,
                Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        }
    }

    KOKKOS_INLINE_FUNCTION
    auto
    getRates(std::size_t reactionId)
    {
        return Kokkos::subview(rates, reactionId, Kokkos::ALL);
    }

    Kokkos::View<double*****, Kokkos::MemoryUnmanaged> productionCoeffs;
    Kokkos::View<double*****, Kokkos::MemoryUnmanaged> dissociationCoeffs;
    Kokkos::View<double**, Kokkos::MemoryUnmanaged> rates;
    ClusterConnectivity<> connectivity;
};
}
}
}
