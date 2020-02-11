#pragma once

#include <Kokkos_View.hpp>
#include <string>

#include <experimental/ReactionNetworkTraits.h>

// TODO: move it
using namespace std::string_literals;

namespace xolotlCore
{
namespace experimental
{
template <typename PlsmContext>
class ClusterCommon;

template <typename TNetwork, typename PlsmContext>
class Cluster;

namespace detail
{
using DefaultMemorySpace = typename Kokkos::View<int*>::traits::memory_space;

template <typename PlsmContext>
struct ContextLabelHelper;

template <>
struct ContextLabelHelper<plsm::OnHost>
{
    static constexpr char label = 'H';
};

template <>
struct ContextLabelHelper<plsm::OnDevice>
{
    static constexpr char label = 'D';
};

template <typename PlsmContext>
constexpr char contextLabel = ContextLabelHelper<PlsmContext>::label;

template <typename TData, typename PlsmContext>
struct ViewTypeHelper;

template <typename TData>
struct ViewTypeHelper<TData, plsm::OnHost>
{
    using DeviceView = Kokkos::View<TData, DefaultMemorySpace>;
    using ViewType = typename DeviceView::HostMirror;
};

template <typename TData>
struct ViewTypeHelper<TData, plsm::OnDevice>
{
    using ViewType = Kokkos::View<TData, DefaultMemorySpace>;
};

template <typename TData, typename PlsmContext>
using ViewType = typename ViewTypeHelper<TData, PlsmContext>::ViewType;

template <typename TView>
struct UnmanagedHelper
{
    using Traits = typename TView::traits;
    using Type =
        Kokkos::View<typename Traits::data_type, typename Traits::array_layout,
            typename Traits::device_type, Kokkos::MemoryUnmanaged>;
};

template <typename TView>
using Unmanaged = typename UnmanagedHelper<TView>::Type;

template <typename PlsmContext = plsm::OnDevice>
struct ClusterDataCommon
{
protected:
    static constexpr char label = contextLabel<PlsmContext>;

public:
    template <typename TData>
    using View = ViewType<TData, PlsmContext>;

    using ClusterType = ClusterCommon<PlsmContext>;

    ClusterDataCommon() = default;

    explicit
    ClusterDataCommon(std::size_t numClusters_, std::size_t gridSize_ = 0)
        :
        numClusters(numClusters_),
        gridSize(gridSize_),
        atomicVolume("Atomic Volume"s + label),
        temperature("Temperature"s + label, gridSize),
        reactionRadius("Reaction Radius"s + label, numClusters),
        formationEnergy("Formation Energy"s + label, numClusters),
        migrationEnergy("Migration Energy"s + label, numClusters),
        diffusionFactor("Diffusion Factor"s + label, numClusters),
        diffusionCoefficient("Diffusion Coefficient"s + label, numClusters,
            gridSize)
    {
    }

    ClusterType
    getCluster(std::size_t clusterId) const noexcept
    {
        return ClusterType(*this, clusterId);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getAtomicVolume() const
    {
        return atomicVolume(0);
    }

    std::size_t numClusters {};
    std::size_t gridSize {};
    View<double[1]> atomicVolume;
    View<double*> temperature;
    View<double*> reactionRadius;
    View<double*> formationEnergy;
    View<double*> migrationEnergy;
    View<double*> diffusionFactor;
    View<double**> diffusionCoefficient;
};

template <typename TNetwork, typename PlsmContext = plsm::OnDevice>
struct ClusterData : ClusterDataCommon<PlsmContext>
{
private:
    using Types = detail::ReactionNetworkTypes<TNetwork>;
    using Props = detail::ReactionNetworkProperties<TNetwork>;
    static constexpr auto nMomentIds = Props::numSpeciesNoI;

public:
    using Superclass = ClusterDataCommon<PlsmContext>;
    using Subpaving = typename Types::Subpaving;
    using TilesView =
        Unmanaged<typename Subpaving::template TilesView<PlsmContext>>;
    using ClusterType = Cluster<TNetwork, PlsmContext>;

    template <typename TData>
    using View = typename Superclass::template View<TData>;

    ClusterData() = default;

    ClusterData(const TilesView& tiles_, std::size_t numClusters_,
            std::size_t gridSize_ = 0)
        :
        Superclass(numClusters_, gridSize_),
        tiles(tiles_),
        momentIds("Moment Ids"s + this->label, numClusters_)
    {
    }

    explicit
    ClusterData(Subpaving& subpaving, std::size_t gridSize_ = 0)
        :
        ClusterData(subpaving.getTiles(PlsmContext{}),
            subpaving.getNumberOfTiles(PlsmContext{}), gridSize_)
    {
    }

    KOKKOS_INLINE_FUNCTION
    ClusterType
    getCluster(std::size_t clusterId) const noexcept
    {
        return ClusterType(*this, clusterId);
    }

    TilesView tiles;
    View<std::size_t*[nMomentIds]> momentIds;
};

template <typename PlsmContext = plsm::OnDevice>
struct ClusterDataCommonRef
{
    using ClusterType = ClusterCommon<PlsmContext>;

    template <typename TData>
    using View = Unmanaged<ViewType<TData, PlsmContext>>;

    ClusterDataCommonRef() = default;

    KOKKOS_INLINE_FUNCTION
    ClusterDataCommonRef(const ClusterDataCommon<PlsmContext>& data)
        :
        numClusters(data.numClusters),
        gridSize(data.gridSize),
        atomicVolume(data.atomicVolume),
        temperature(data.temperature),
        reactionRadius(data.reactionRadius),
        formationEnergy(data.formationEnergy),
        migrationEnergy(data.migrationEnergy),
        diffusionFactor(data.diffusionFactor),
        diffusionCoefficient(data.diffusionCoefficient)
    {
    }

    KOKKOS_INLINE_FUNCTION
    ClusterType
    getCluster(std::size_t clusterId) const noexcept
    {
        return ClusterType(*this, clusterId);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getAtomicVolume() const
    {
        return atomicVolume(0);
    }

    std::size_t numClusters {};
    std::size_t gridSize {};
    View<double[1]> atomicVolume;
    View<double*> temperature;
    View<double*> reactionRadius;
    View<double**> diffusionCoefficient;
    View<double*> formationEnergy;
    View<double*> migrationEnergy;
    View<double*> diffusionFactor;
};

template <typename TNetwork, typename PlsmContext = plsm::OnDevice>
struct ClusterDataRef : ClusterDataCommonRef<PlsmContext>
{
private:
    using Types = detail::ReactionNetworkTypes<TNetwork>;
    using Props = detail::ReactionNetworkProperties<TNetwork>;
    static constexpr auto nMomentIds = Props::numSpeciesNoI;

public:
    using Superclass = ClusterDataCommonRef<PlsmContext>;
    using Subpaving = typename Types::Subpaving;
    using TilesView =
        Unmanaged<typename Subpaving::template TilesView<PlsmContext>>;
    using ClusterType = Cluster<TNetwork, PlsmContext>;

    template <typename TData>
    using View = typename Superclass::template View<TData>;

    ClusterDataRef() = default;

    KOKKOS_INLINE_FUNCTION
    ClusterDataRef(const ClusterData<TNetwork, PlsmContext>& data)
        :
        Superclass(data),
        tiles(data.tiles),
        momentIds(data.momentIds)
    {
    }

    KOKKOS_INLINE_FUNCTION
    ClusterType
    getCluster(std::size_t clusterId) const noexcept
    {
        return ClusterType(*this, clusterId);
    }

    TilesView tiles;
    View<std::size_t*[nMomentIds]> momentIds;
};
}
}
}
