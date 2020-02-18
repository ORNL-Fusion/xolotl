#pragma once

#include <Kokkos_Core.hpp>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
struct ReactionData
{
    ReactionData() = default;

    ReactionData(std::size_t numProductionReactions,
            std::size_t numDissociationReactions, std::size_t numDOFs,
            std::size_t numSpeciesNoI, std::size_t gridSize)
        :
        coeffExtent(numSpeciesNoI + 1),
        numReactions(numProductionReactions + numDissociationReactions),
        productionCoeffs("Production Coefficients",
            numProductionReactions, coeffExtent, coeffExtent, 4, coeffExtent),
        dissociationCoeffs("Dissociation Coefficients",
            numDissociationReactions, coeffExtent, 1, 3, coeffExtent),
        rates("Reaction Rates", numReactions, gridSize),
        inverseMap(
            Kokkos::ViewAllocateWithoutInitializing("Reactions Inverse Map"),
            numDOFs, numDOFs)
    {
        initializeInverseMap(numDOFs);
    }

    void
    initializeInverseMap(std::size_t numDOFs)
    {
        auto invMap = inverseMap;
        using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        Kokkos::parallel_for(Range2D({0, 0}, {numDOFs, numDOFs}),
                KOKKOS_LAMBDA (std::size_t i, std::size_t j) {
            invMap(i,j) = plsm::invalid<std::size_t>;
        });
    }

    std::size_t coeffExtent {};
    std::size_t numReactions {};
    Kokkos::View<double*****> productionCoeffs;
    Kokkos::View<double*****> dissociationCoeffs;
    Kokkos::View<double**> rates;

    // TODO: the original code uses an actual map here because it is sparse
    //       Does this need to be Kokkos::Crs?
    Kokkos::View<std::size_t**> inverseMap;
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
    Kokkos::View<std::size_t**, Kokkos::MemoryUnmanaged> inverseMap;
};
}
}
}
