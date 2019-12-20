#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TSubpaving, typename TMemSpace>
class Cluster;

namespace detail
{
//TODO: Change plsm:: OnHost/OnDevice to be aliases for Kokkos::HostSpace and
//this
using DefaultMemorySpace = typename Kokkos::View<int*>::traits::memory_space;

template <typename PLSMContext>
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

template <typename PLSMContext>
using KokkosMemorySpace = typename KokkosMemorySpaceHelper<PLSMContext>::Space;

template <typename TMemSpace>
struct ContextHelper
{
    using Context = plsm::OnDevice;
};

template <>
struct ContextHelper<Kokkos::HostSpace>
{
    using Context = plsm::OnHost;
};

template <typename TMemSpace>
using Context = typename ContextHelper<TMemSpace>::Context;

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

template <typename TSubpaving, typename TMemSpace = DefaultMemorySpace>
struct ClusterData
{
    using MemSpace = TMemSpace;
    using Subpaving = TSubpaving;
    //TODO: Apply MemSpace
    using TilesView =
        Unmanaged<typename Subpaving::template TilesView<plsm::OnDevice>>;
    using ClusterType = Cluster<Subpaving, detail::Context<MemSpace>>;

    template <typename TData>
    using View = Kokkos::View<TData, MemSpace>;

    ClusterData() = default;

    explicit
    ClusterData(std::size_t numClusters_, std::size_t gridSize_ = 0)
        :
        numClusters(numClusters_),
        gridSize(gridSize_),
        temperature("Temperature", gridSize),
        momentIds("Moment Ids", numClusters),
        reactionRadius("Reaction Radius", numClusters),
        formationEnergy("Formation Energy", numClusters),
        migrationEnergy("Migration Energy", numClusters),
        diffusionFactor("Diffusion Factor", numClusters),
        diffusionCoefficient("Diffusion Coefficient", numClusters, gridSize)
    {
    }

    ClusterData(const ClusterData&) = default;
    ClusterData& operator=(const ClusterData&) = default;

    KOKKOS_INLINE_FUNCTION
    ClusterType
    getCluster(std::size_t clusterId)
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

template <typename TSubpaving, typename TMemSpace = DefaultMemorySpace>
struct ClusterDataRef
{
    using MemSpace = TMemSpace;
    using Subpaving = TSubpaving;

    //TODO: Apply MemSpace
    using TilesView =
        Unmanaged<typename Subpaving::template TilesView<plsm::OnDevice>>;

    using ClusterType = Cluster<Subpaving, detail::Context<MemSpace>>;

    template <typename TData>
    using View = Kokkos::View<TData, MemSpace, Kokkos::MemoryUnmanaged>;

    KOKKOS_INLINE_FUNCTION
    ClusterDataRef() = default;

    KOKKOS_INLINE_FUNCTION
    ClusterDataRef(detail::ClusterData<Subpaving, MemSpace>& data)
        :
        numClusters(data.numClusters),
        gridSize(data.gridSize),
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
    ClusterDataRef(const ClusterDataRef&) = default;
    KOKKOS_INLINE_FUNCTION
    ClusterDataRef& operator=(const ClusterDataRef&) = default;

    KOKKOS_INLINE_FUNCTION
    ClusterType
    getCluster(std::size_t clusterId)
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


template <typename TSubpaving, typename PLSMContext>
class Cluster
{
public:
    using MemSpace = detail::KokkosMemorySpace<PLSMContext>;
    using Subpaving = TSubpaving;
    using Region = typename Subpaving::RegionType;
    using ClusterData = detail::ClusterData<Subpaving, MemSpace>;
    using ClusterDataRef = detail::ClusterDataRef<Subpaving, MemSpace>;

    Cluster() = delete;

    KOKKOS_INLINE_FUNCTION
    Cluster(ClusterData& data, std::size_t id)
        :
        _data{data},
        _id{id}
    {
    }

    KOKKOS_INLINE_FUNCTION
    Cluster(ClusterDataRef& data, std::size_t id)
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
