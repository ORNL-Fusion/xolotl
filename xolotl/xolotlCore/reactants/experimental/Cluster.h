#pragma once

#include <typeinfo>

namespace xolotlCore
{
namespace experimental
{
template <typename TSubpaving, typename TMemSpace>
class Cluster;

namespace detail
{
using DefaultMemorySpace = typename Kokkos::View<int*>::traits::memory_space;

template <typename PlsmContext>
struct KokkosMemorySpaceHelper;

template <>
struct KokkosMemorySpaceHelper<plsm::OnHost>
{
    using Space = Kokkos::HostSpace;
};

template <>
struct KokkosMemorySpaceHelper<plsm::OnDevice>
{
    using Space = DefaultMemorySpace;
};

template <typename PlsmContext>
using KokkosMemorySpace = typename KokkosMemorySpaceHelper<PlsmContext>::Space;

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

template <typename TSubpaving, typename PlsmContext = plsm::OnDevice>
struct ClusterData
{
private:
    static constexpr char label = contextLabel<PlsmContext>;

public:
    using Subpaving = TSubpaving;
    using TilesView =
        Unmanaged<typename Subpaving::template TilesView<PlsmContext>>;
    using ClusterType = Cluster<Subpaving, PlsmContext>;

    template <typename TData>
    using View = ViewType<TData, PlsmContext>;

    ClusterData() = default;

    ClusterData(const TilesView& tiles_, std::size_t numClusters_,
            std::size_t gridSize_ = 0)
        :
        numClusters(numClusters_),
        gridSize(gridSize_),
        atomicVolume("Atomic Volume"s + label),
        temperature("Temperature"s + label, gridSize),
        tiles(tiles_),
        momentIds("Moment Ids"s + label, numClusters),
        reactionRadius("Reaction Radius"s + label, numClusters),
        formationEnergy("Formation Energy"s + label, numClusters),
        migrationEnergy("Migration Energy"s + label, numClusters),
        diffusionFactor("Diffusion Factor"s + label, numClusters),
        diffusionCoefficient("Diffusion Coefficient"s + label, numClusters,
            gridSize)
    {
    }

    explicit
    ClusterData(Subpaving& subpaving, std::size_t gridSize_ = 0)
        :
        ClusterData(subpaving.getTiles(PlsmContext{}),
            subpaving.getNumberOfTiles(PlsmContext{}), gridSize_)
    {
    }

    ClusterData(const ClusterData&) = default;
    ClusterData& operator=(const ClusterData&) = default;

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
    TilesView tiles;
    View<std::size_t*[4]> momentIds;
    View<double*> reactionRadius;
    View<double*> formationEnergy;
    View<double*> migrationEnergy;
    View<double*> diffusionFactor;
    View<double**> diffusionCoefficient;
};

template <typename TSubpaving, typename PlsmContext = plsm::OnDevice>
struct ClusterDataRef
{
    using Subpaving = TSubpaving;

    using TilesView =
        Unmanaged<typename Subpaving::template TilesView<PlsmContext>>;

    using ClusterType = Cluster<Subpaving, PlsmContext>;

    template <typename TData>
    using View = Unmanaged<ViewType<TData, PlsmContext>>;

    KOKKOS_INLINE_FUNCTION
    ClusterDataRef() = default;

    KOKKOS_INLINE_FUNCTION
    ClusterDataRef(const detail::ClusterData<Subpaving, PlsmContext>& data)
        :
        numClusters(data.numClusters),
        gridSize(data.gridSize),
        atomicVolume(data.atomicVolume),
        temperature(data.temperature),
        tiles(data.tiles),
        momentIds(data.momentIds),
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
    TilesView tiles;
    View<std::size_t*[4]> momentIds;
    View<double*> reactionRadius;
    View<double**> diffusionCoefficient;
    View<double*> formationEnergy;
    View<double*> migrationEnergy;
    View<double*> diffusionFactor;
};
}


template <typename TSubpaving, typename PlsmContext>
class Cluster
{
public:
    using Subpaving = TSubpaving;
    using Region = typename Subpaving::RegionType;
    using ClusterData = detail::ClusterData<Subpaving, PlsmContext>;
    using ClusterDataRef = detail::ClusterDataRef<Subpaving, PlsmContext>;

    Cluster() = delete;

    KOKKOS_INLINE_FUNCTION
    Cluster(const ClusterData& data, std::size_t id)
        :
        _data{data},
        _id{id}
    {
    }

    KOKKOS_INLINE_FUNCTION
    Cluster(const ClusterDataRef& data, std::size_t id)
        :
        _data{data},
        _id{id}
    {
    }

    KOKKOS_INLINE_FUNCTION
    Region
    getRegion() const
    {
        return _data.tiles(_id).getRegion();
    }

    KOKKOS_INLINE_FUNCTION
    decltype(auto)
    getMomentIds()
    {
        return Kokkos::subview(_data.momentIds, _id, Kokkos::ALL);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getReactionRadius()
    {
        return _data.reactionRadius(_id);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getFormationEnergy()
    {
        return _data.formationEnergy(_id);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getDiffusionCoefficient(std::size_t gridIndex)
    {
        return _data.diffusionCoefficient(_id, gridIndex);
    }

    std::size_t
    getId()
    {
        return _id;
    }

private:
    ClusterDataRef _data;
    std::size_t _id {plsm::invalid<std::size_t>};
};
}
}
