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
struct ReactionInverseMap
{
    using Connectivity = Kokkos::Crs<std::size_t, TMemSpace>;
    using HostMirror =
        ReactionInverseMap<
            typename Kokkos::View<int*, TMemSpace>::traits::host_mirror_space>;

    Connectivity connectivity;

    KOKKOS_INLINE_FUNCTION
    std::size_t
    operator()(std::size_t rowId, std::size_t columnId) const
    {
        for (auto pos = connectivity.row_map(rowId);
                pos < connectivity.row_map(rowId+1); ++pos) {
            if (connectivity.entries(pos) == columnId) {
                return pos;
            }
        }
        return plsm::invalid<std::size_t>;
    }
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

    std::size_t coeffExtent {};
    std::size_t numReactions {};
    Kokkos::View<double*****> productionCoeffs;
    Kokkos::View<double*****> dissociationCoeffs;
    Kokkos::View<double**> rates;
    ReactionInverseMap<> inverseMap;
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
        inverseMap(data.inverseMap)
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
    ReactionInverseMap<> inverseMap;
};
}
}
}
