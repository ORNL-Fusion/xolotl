#pragma once

#include <type_traits>

#include <Kokkos_Core.hpp>

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#define DEVICE_LAMBDA [=] __device__
#else
#define DEVICE_LAMBDA [=]
#endif

#include <xolotl/core/reactants/ReactionData.h>
#include <xolotl/core/reactants/ReactionNetworkTraits.h>
#include <xolotl/core/reactants/TupleUtility.h>

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

    void
    allocateCoefficients(CoefficientsView& coeffs) const
    {
        coeffs = ReactionType::allocateCoefficientsView(_numReactions);
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

template <std::size_t NumReactionTypes, typename... TReactions>
class ReactionSetMixinChain
{
    template <std::size_t, typename...>
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
    allocateCoefficients(Kokkos::Array<CoefficientsView, NumReactionTypes>&)
        const
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

    void
    getReactionBeginIndices(Kokkos::Array<IndexType, NumReactionTypes+1>& ids)
        const
    {
        ids[NumReactionTypes] = _indexBegin;
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
    KOKKOS_INLINE_FUNCTION
    IndexType
    getIndexBegin() const noexcept
    {
        return _indexBegin;
    }

private:
    IndexType _indexBegin {};
};

template <std::size_t NumReactionTypes, typename TReaction,
    typename... TOtherReactions>
class ReactionSetMixinChain<NumReactionTypes, TReaction, TOtherReactions...> :
    public
        ReactionSet<TReaction>,
        ReactionSetMixinChain<NumReactionTypes, TOtherReactions...>
{
    template <std::size_t, typename...>
    friend class ReactionSetMixinChain;

public:
    using IndexType = ReactionNetworkIndexType;
    using Head = ReactionSet<TReaction>;
    using Tail = ReactionSetMixinChain<NumReactionTypes, TOtherReactions...>;

private:
    static constexpr std::size_t
    getReactionTypeIndex()
    {
        return NumReactionTypes - (sizeof...(TOtherReactions) + 1);
    }

public:
    ReactionSetMixinChain() = default;

    template <typename THeadView, typename... TTailViews>
    ReactionSetMixinChain(IndexType indexBegin, THeadView view,
            TTailViews... views)
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

    void
    allocateCoefficients(
        Kokkos::Array<CoefficientsView, NumReactionTypes>& coeffs) const
    {
        Head::allocateCoefficients(coeffs[getReactionTypeIndex()]);
        Tail::allocateCoefficients(coeffs);
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

    IndexType
    getNumberOfRates() const noexcept
    {
        return Head::getNumberOfRates() + Tail::getNumberOfRates();
    }

    void
    getReactionBeginIndices(Kokkos::Array<IndexType, NumReactionTypes+1>& ids)
        const
    {
        ids[getReactionTypeIndex()] = _indexBegin;
        Tail::getReactionBeginIndices(ids);
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
    KOKKOS_INLINE_FUNCTION
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
    using Type = ReactionSetMixinChain<sizeof...(TReactions), TReactions...>;
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

private:
    static constexpr std::size_t numReactionTypes =
        std::tuple_size<ReactionTypes>::value;

    static constexpr std::size_t numSpeciesNoI =
        NetworkType::getNumberOfSpeciesNoI();

public:
    ReactionCollection() = default;

    template <typename... TViews>
    ReactionCollection(IndexType gridSize, TViews... views)
        :
        _chain(views...),
        _data(_chain.getNumberOfReactions(), gridSize,
            getReactionBeginIndices())
    {
        static_assert(sizeof...(TViews) == numReactionTypes,
            "Construction from views requires the number of views to match the "
            "number of reaction types in the ReactionCollection");

        _chain.allocateCoefficients(_data.coeffs);
    }

    void
    setGridSize(IndexType gridSize)
    {
        _data.setGridSize(gridSize);
    }

    const ClusterConnectivity<>&
    getConnectivity() const
    {
        return _data.connectivity;
    }

    void
    setConnectivity(const ClusterConnectivity<>& connectivity)
    {
        _data.connectivity = connectivity;
    }

    std::uint64_t
    getDeviceMemorySize() const noexcept
    {
        std::uint64_t ret = _chain.getDeviceMemorySize();
        ret += _data.getDeviceMemorySize();
        return ret;
    }

    IndexType
    getNumberOfReactions() const noexcept
    {
        return _data.numReactions;
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
        _data.numReactions = _chain.getNumberOfReactions();
    }

    void
    constructAll(ClusterDataRef clusterData,
        Kokkos::View<ClusterSet*> clusterSets)
    {
        auto chain = _chain;
        auto reactionData = ReactionDataRef<NetworkType>(_data);
        Kokkos::parallel_for(_data.numReactions,
                DEVICE_LAMBDA (const IndexType i) {
            chain.apply(DEVICE_LAMBDA (auto& reaction) {
                using ReactionType =
                    std::remove_reference_t<decltype(reaction)>;
                reaction =
                    ReactionType(reactionData, clusterData, i, clusterSets(i));
            }, i);
        });
    }

    void
    updateAll(ClusterDataRef clusterData)
    {
        auto reactionData = ReactionDataRef<NetworkType>(_data);
        apply(DEVICE_LAMBDA (auto&& reaction) {
            reaction.updateData(reactionData, clusterData);
        });
    }

    double
    getLargestRate() const
    {
        using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        double largestRate = 0.0;
        auto rates = _data.rates;
        auto nRates = rates.extent(0);
        auto gridSize = rates.extent(1);
        Kokkos::parallel_reduce(Range2D({0, 0}, {nRates, gridSize}),
                KOKKOS_LAMBDA (const IndexType i, const IndexType j,
                    double& max) {
            if (rates(i, j) > max) {
                max = rates(i, j);
            }
        }, Kokkos::Max<double>(largestRate));
        Kokkos::fence();
        return largestRate;
    }

    template <typename F>
    void
    apply(const F& func)
    {
        auto chain = _chain;
        Kokkos::parallel_for(_data.numReactions, DEVICE_LAMBDA (const IndexType i) {
            chain.apply(func, i);
        });
    }

    template <typename F, typename T>
    void
    reduce(const F& func, T& out)
    {
        auto chain = _chain;
        Kokkos::parallel_reduce(_data.numReactions,
                DEVICE_LAMBDA (const IndexType i, T& local) {
            chain.reduce(func, i, local);
        }, out);
    }

private:
    Kokkos::Array<IndexType, numReactionTypes+1>
    getReactionBeginIndices() const
    {
        Kokkos::Array<IndexType, numReactionTypes+1> ret;
        _chain.getReactionBeginIndices(ret);
        return ret;
    }

private:
    ReactionSetChain<ReactionTypes> _chain;
    ReactionData<NetworkType> _data;
};
}
}
}
