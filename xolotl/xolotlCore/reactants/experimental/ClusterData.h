#pragma once

#include <string>

#include <Kokkos_View.hpp>

#include <experimental/ReactionNetworkTraits.h>
#include <experimental/MemorySpace.h>

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

inline
std::string
labelStr(const char labelChar)
{
    return std::string("_") + labelChar;
}

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
    using IndexType = detail::ReactionNetworkIndexType;

    ClusterDataCommon() = default;

    explicit
    ClusterDataCommon(IndexType numClusters_, IndexType gridSize_ = 0)
        :
        numClusters(numClusters_),
        gridSize(gridSize_),
        atomicVolume("Atomic Volume" + labelStr(label)),
        latticeParameter("Lattice Parameter" + labelStr(label)),
        fissionRate("Fission Rate" + labelStr(label)),
        zeta("Zeta" + labelStr(label)),
        enableStdReaction("Enable Std Reaction" + labelStr(label)),
        enableReSolution("Enable Re-Solution Process" + labelStr(label)),
        temperature("Temperature" + labelStr(label), gridSize),
        reactionRadius("Reaction Radius" + labelStr(label), numClusters),
        formationEnergy("Formation Energy" + labelStr(label), numClusters),
        migrationEnergy("Migration Energy" + labelStr(label), numClusters),
        diffusionFactor("Diffusion Factor" + labelStr(label), numClusters),
        diffusionCoefficient("Diffusion Coefficient" + labelStr(label),
            numClusters, gridSize)
    {
    }

    ClusterType
    getCluster(IndexType clusterId) const noexcept
    {
        return ClusterType(*this, clusterId);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getAtomicVolume() const
    {
        return atomicVolume(0);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getLatticeParamter() const
    {
        return latticeParameter(0);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getFissionRate() const
    {
        return fissionRate(0);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getZeta() const
    {
        return zeta(0);
    }

    KOKKOS_INLINE_FUNCTION
    bool
    getEnableStdReaction() const
    {
        return enableStdReaction(0);
    }

    KOKKOS_INLINE_FUNCTION
    bool
    getEnableReSolution() const
    {
        return enableReSolution(0);
    }

    void
    setGridSize(IndexType gridSize_) {
        gridSize = gridSize_;
        temperature = View<double*>("Temperature" + labelStr(label), gridSize);
        diffusionCoefficient = View<double**>("Diffusion Coefficient" + labelStr(label),
                numClusters, gridSize);
    }

    IndexType numClusters {};
    IndexType gridSize {};
    View<double[1]> atomicVolume;
    View<double[1]> latticeParameter;
    View<double[1]> fissionRate;
    View<double[1]> zeta;
    View<bool[1]> enableStdReaction;
    View<bool[1]> enableReSolution;
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
    using IndexType = typename Types::IndexType;

    template <typename TData>
    using View = typename Superclass::template View<TData>;

    ClusterData() = default;

    ClusterData(const TilesView& tiles_, IndexType numClusters_,
            IndexType gridSize_ = 0)
        :
        Superclass(numClusters_, gridSize_),
        tiles(tiles_),
        momentIds("Moment Ids" + labelStr(this->label), numClusters_)
    {
    }

    explicit
    ClusterData(Subpaving& subpaving, IndexType gridSize_ = 0)
        :
        ClusterData(subpaving.getTiles(PlsmContext{}),
            subpaving.getNumberOfTiles(PlsmContext{}), gridSize_)
    {
    }

    KOKKOS_INLINE_FUNCTION
    ClusterType
    getCluster(IndexType clusterId) const noexcept
    {
        return ClusterType(*this, clusterId);
    }

    TilesView tiles;
    View<IndexType*[nMomentIds]> momentIds;
};

template <typename PlsmContext = plsm::OnDevice>
struct ClusterDataCommonRef
{
    using ClusterType = ClusterCommon<PlsmContext>;
    using IndexType = detail::ReactionNetworkIndexType;

    template <typename TData>
    using View = Unmanaged<ViewType<TData, PlsmContext>>;

    ClusterDataCommonRef() = default;

    KOKKOS_INLINE_FUNCTION
    ClusterDataCommonRef(const ClusterDataCommon<PlsmContext>& data)
        :
        numClusters(data.numClusters),
        gridSize(data.gridSize),
        atomicVolume(data.atomicVolume),
        latticeParameter(data.latticeParameter),
        fissionRate(data.fissionRate),
        zeta(data.zeta),
        enableStdReaction(data.enableStdReaction),
        enableReSolution(data.enableReSolution),
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
    getCluster(IndexType clusterId) const noexcept
    {
        return ClusterType(*this, clusterId);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getAtomicVolume() const
    {
        return atomicVolume(0);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getLatticeParameter() const
    {
        return latticeParameter(0);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getFissionRate() const
    {
        return fissionRate(0);
    }

    KOKKOS_INLINE_FUNCTION
    double
    getZeta() const
    {
        return zeta(0);
    }

    KOKKOS_INLINE_FUNCTION
    bool
    getEnableStdReaction() const
    {
        return enableStdReaction(0);
    }

    KOKKOS_INLINE_FUNCTION
    bool
    getEnableReSolution() const
    {
        return enableReSolution(0);
    }

    IndexType numClusters {};
    IndexType gridSize {};
    View<double[1]> atomicVolume;
    View<double[1]> latticeParameter;
    View<double[1]> fissionRate;
    View<double[1]> zeta;
    View<bool[1]> enableStdReaction;
    View<bool[1]> enableReSolution;
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
    using IndexType = typename Types::IndexType;

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
    getCluster(IndexType clusterId) const noexcept
    {
        return ClusterType(*this, clusterId);
    }

    TilesView tiles;
    View<IndexType*[nMomentIds]> momentIds;
};
}
}
}
