#pragma once

#include <Kokkos_Core.hpp>

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#define DEVICE_LAMBDA [=] __device__
#else
#define DEVICE_LAMBDA [=]
#endif

#include <experimental/ReactionNetworkTraits.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TTraits>
class ReactionCollection
{
public:
    using Traits = TTraits;
    using IndexType = detail::ReactionNetworkIndexType;
    using ProductionReactionType = typename Traits::ProductionReactionType;
    using DissociationReactionType = typename Traits::DissociationReactionType;

    ReactionCollection() = default;

    ReactionCollection(Kokkos::View<ProductionReactionType*> prodReactions,
            Kokkos::View<DissociationReactionType*> dissReactions)
        :
        _prodReactions(prodReactions),
        _dissReactions(dissReactions),
        _numProdReactions(prodReactions.size()),
        _numReactions(_numProdReactions + dissReactions.size())
    {
    }

    std::uint64_t
    getDeviceMemorySize() const noexcept
    {
        std::uint64_t ret {};
        ret += _prodReactions.required_allocation_size(_prodReactions.size());
        ret += _dissReactions.required_allocation_size(_dissReactions.size());
        return ret;
    }

    template <typename F>
    void
    apply(const F& func)
    {
        auto prodReactions = _prodReactions;
        auto dissReactions = _dissReactions;
        auto numProdReactions = _numProdReactions;
        Kokkos::parallel_for(_numReactions, DEVICE_LAMBDA (const IndexType i) {
            if (i < numProdReactions) {
                func(prodReactions(i));
            }
            else {
                func(dissReactions(i - numProdReactions));
            }
        });
    }

    template <typename F, typename T>
    void
    reduce(const F& func, T& out)
    {
        auto prodReactions = _prodReactions;
        auto dissReactions = _dissReactions;
        auto numProdReactions = _numProdReactions;
        Kokkos::parallel_reduce(_numReactions,
                DEVICE_LAMBDA (const IndexType i, T& local) {
            if (i < numProdReactions) {
                func(prodReactions(i), local);
            }
            else {
                func(dissReactions(i - numProdReactions), local);
            }
        }, out);
    }

private:
    Kokkos::View<ProductionReactionType*> _prodReactions;
    Kokkos::View<DissociationReactionType*> _dissReactions;
    IndexType _numProdReactions;
    IndexType _numReactions;
};
}
}
}
