#pragma once

#include <type_traits>

#include <Kokkos_Core.hpp>

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#define DEVICE_LAMBDA [=] __device__
#else
#define DEVICE_LAMBDA [=]
#endif

#include <experimental/ReactionNetworkTraits.h>
#include <experimental/TupleUtility.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TReaction>
class ReactionSet
{
public:
    using ReactionType = TReaction;
    using IndexType = ReactionNetworkIndexType;

    ReactionSet() = default;

    ReactionSet(Kokkos::View<ReactionType*> rView)
        :
        _reactions(rView),
        _numReactions(rView.size())
    {
    }

    std::uint64_t
    getDeviceMemorySize() const noexcept
    {
        return _reactions.required_allocation_size(_numReactions);
    }

    IndexType
    getNumberOfReactions() const noexcept
    {
        return _numReactions;
    }

    Kokkos::View<ReactionType*>
    getView() const
    {
        return _reactions;
    }

    void
    setView(Kokkos::View<ReactionType*> view)
    {
        _reactions = view;
        _numReactions = view.size();
    }

    template <typename F>
    KOKKOS_INLINE_FUNCTION
    void
    apply(const F& func, const IndexType i) const
    {
        func(_reactions(i));
    }

    template <typename F, typename T>
    KOKKOS_INLINE_FUNCTION
    void
    reduce(const F& func, const IndexType i, T& local) const
    {
        func(_reactions(i), local);
    }

private:
    Kokkos::View<ReactionType*> _reactions;
    IndexType _numReactions {};
};

template <typename... TReactions>
class ReactionSetMixinChain
{
    template <typename...>
    friend class ReactionSetMixinChain;

public:
    using IndexType = ReactionNetworkIndexType;

    ReactionSetMixinChain() = default;

    ReactionSetMixinChain(IndexType indexBegin)
        :
        _indexBegin(indexBegin)
    {
    }

    void
    updateIndices(IndexType)
    {
    }

    void
    getView() const
    {
    }

    std::uint64_t
    getDeviceMemorySize() const noexcept
    {
        return 0;
    }

    ReactionNetworkIndexType
    getNumberOfReactions() const noexcept
    {
        return 0;
    }

    template <typename F>
    KOKKOS_INLINE_FUNCTION
    void
    apply(const F&, const IndexType) const noexcept
    {
    }

    template <typename F, typename T>
    KOKKOS_INLINE_FUNCTION
    void
    reduce(const F&, const IndexType, T&) const
    {
    }

private:
    IndexType
    getIndexBegin() const noexcept
    {
        return _indexBegin;
    }

private:
    IndexType _indexBegin {};
};

template <typename TReaction, typename... TOtherReactions>
class ReactionSetMixinChain<TReaction, TOtherReactions...> :
    public ReactionSet<TReaction>, ReactionSetMixinChain<TOtherReactions...>
{
    template <typename...>
    friend class ReactionSetMixinChain;

public:
    using IndexType = ReactionNetworkIndexType;
    using Head = ReactionSet<TReaction>;
    using Tail = ReactionSetMixinChain<TOtherReactions...>;

    ReactionSetMixinChain() = default;

    template <typename THeadView, typename... TTailViews>
    ReactionSetMixinChain(IndexType indexBegin, THeadView view, TTailViews... views)
        :
        Head(view),
        Tail(indexBegin + view.size(), views...),
        _indexBegin(indexBegin)
    {
    }

    template <typename... TViews>
    ReactionSetMixinChain(TViews... views)
        :
        ReactionSetMixinChain(static_cast<IndexType>(0), views...)
    {
    }

    void
    updateIndices(IndexType indexBegin)
    {
        _indexBegin = indexBegin;
        Tail::updateIndices(indexBegin + Head::getNumberOfReactions());
    }

    template <typename TGetReaction,
        std::enable_if_t<std::is_same<TGetReaction, TReaction>::value,
            int> = 0>
    decltype(auto)
    getView() const
    {
        return Head::getView();
    }

    template <typename TGetReaction,
        std::enable_if_t<!std::is_same<TGetReaction, TReaction>::value,
            int> = 0>
    decltype(auto)
    getView() const
    {
        return Tail::template getView<TGetReaction>();
    }

    template <typename TSetReaction,
        std::enable_if_t<std::is_same<TSetReaction, TReaction>::value,
            int> = 0>
    void
    setView(Kokkos::View<TSetReaction*> view)
    {
        Head::setView(view);
        updateIndices(_indexBegin);
    }

    template <typename TSetReaction,
        std::enable_if_t<!std::is_same<TSetReaction, TReaction>::value,
            int> = 0>
    void
    setView(Kokkos::View<TSetReaction*> view)
    {
        Tail::setView(view);
    }

    std::uint64_t
    getDeviceMemorySize() const noexcept
    {
        return Head::getDeviceMemorySize() + Tail::getDeviceMemorySize();
    }

    IndexType
    getNumberOfReactions() const noexcept
    {
        return Head::getNumberOfReactions() + Tail::getNumberOfReactions();
    }

    template <typename F>
    KOKKOS_INLINE_FUNCTION
    void
    apply(const F& func, const IndexType i) const
    {
        if (i < Tail::getIndexBegin()) {
            Head::apply(func, i - _indexBegin);
        }
        else {
            Tail::apply(func, i);
        }
    }

    template <typename F, typename T>
    KOKKOS_INLINE_FUNCTION
    void
    reduce(const F& func, const IndexType i, T& local) const
    {
        if (i < Tail::getIndexBegin()) {
            Head::reduce(func, i - _indexBegin, local);
        }
        else {
            Tail::reduce(func, i, local);
        }
    }

private:
    IndexType
    getIndexBegin() const noexcept
    {
        return _indexBegin;
    }

private:
    IndexType _indexBegin {};
};

template <typename TReactionTypeList>
struct ReactionSetChainHelper;

template <typename... TReactions>
struct ReactionSetChainHelper<std::tuple<TReactions...>>
{
    using Type = ReactionSetMixinChain<TReactions...>;
};

template <typename TReactionTypeList>
using ReactionSetChain =
    typename ReactionSetChainHelper<TReactionTypeList>::Type;

template <typename TNetwork>
class ReactionCollection
{
public:
    using NetworkType = TNetwork;
    using IndexType = detail::ReactionNetworkIndexType;
    using ReactionTypes = ReactionTypeList<NetworkType>;
    using Types = ReactionNetworkTypes<NetworkType>;
    using ClusterDataRef = typename Types::ClusterDataRef;

    ReactionCollection() = default;

    template <typename... TViews>
    ReactionCollection(TViews... views)
        :
        _chain(views...),
        _numReactions(_chain.getNumberOfReactions())
    {
        static_assert(
            sizeof...(TViews) == std::tuple_size<ReactionTypes>::value,
            "Construction from views requires the number of views to match the "
            "number of reaction types in the ReactionCollection");
    }

    std::uint64_t
    getDeviceMemorySize() const noexcept
    {
        return _chain.getDeviceMemorySize();
    }

    IndexType
    getNumberOfReactions() const noexcept
    {
        return _numReactions;
    }

    template <typename TReaction>
    Kokkos::View<TReaction*>
    getView() const
    {
        return _chain.template getView<TReaction>();
    }

    template <typename TReaction>
    void
    setView(Kokkos::View<TReaction*> view)
    {
        _chain.setView(view);
        _numReactions = _chain.getNumberOfReactions();
    }

    void
    construct(ReactionDataRef reactionData, ClusterDataRef clusterData,
        Kokkos::View<ClusterSet*> clusterSets)
    {
        auto chain = _chain;
        Kokkos::parallel_for(_numReactions, DEVICE_LAMBDA (const IndexType i) {
            chain.apply(DEVICE_LAMBDA (auto& reaction) {
                using ReactionType =
                    std::remove_reference_t<decltype(reaction)>;
                reaction =
                    ReactionType(reactionData, clusterData, i, clusterSets(i));
            }, i);
        });
    }

    template <typename F>
    void
    apply(const F& func)
    {
        auto chain = _chain;
        Kokkos::parallel_for(_numReactions, DEVICE_LAMBDA (const IndexType i) {
            chain.apply(func, i);
        });
    }

    template <typename F, typename T>
    void
    reduce(const F& func, T& out)
    {
        auto chain = _chain;
        Kokkos::parallel_reduce(_numReactions,
                DEVICE_LAMBDA (const IndexType i, T& local) {
            chain.reduce(func, i, local);
        }, out);
    }

private:
    ReactionSetChain<ReactionTypes> _chain;
    IndexType _numReactions {};
};
}
}
}
