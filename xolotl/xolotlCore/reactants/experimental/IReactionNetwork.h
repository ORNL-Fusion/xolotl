#pragma once

#include <unordered_map>

#include <Kokkos_Core.hpp>
#include <Kokkos_Crs.hpp>

#include <experimental/ClusterData.h>
#include <experimental/Cluster.h>

namespace xolotlCore
{
namespace experimental
{

class IReactionNetwork
{
public:
    using AmountType = detail::ReactionNetworkAmountType;
    using ConcentrationsView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
    using OwnedConcentrationsView = Kokkos::View<double*>;
    using FluxesView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
    using OwnedFluxesView = Kokkos::View<double*>;
    using Connectivity = Kokkos::Crs<std::size_t, detail::DefaultMemorySpace>;
    using SparseFillMap = std::unordered_map<int, std::vector<int>>;

    IReactionNetwork(std::size_t gridSize)
        :
        _gridSize(gridSize)
    {
    }

    virtual
    ~IReactionNetwork()
    {
    }

    virtual std::uint64_t
    getDeviceMemorySize() const noexcept
    {
        return 0;
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t
    getDOF() const noexcept
    {
        return _numDOFs;
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t
    getNumClusters() const noexcept
    {
        return _numClusters;
    }

    KOKKOS_INLINE_FUNCTION
    double
    getLatticeParameter() const noexcept
    {
        return _latticeParameter;
    }

    virtual void
    setLatticeParameter(double latticeParameter) = 0;

    KOKKOS_INLINE_FUNCTION
    double
    getAtomicVolume() const noexcept
    {
        return _atomicVolume;
    }

    KOKKOS_INLINE_FUNCTION
    double
    getInterstitialBias() const noexcept
    {
        return _interstitialBias;
    }

    void
    setInterstitialBias(double interstitialBias) noexcept
    {
        _interstitialBias = interstitialBias;
    }

    KOKKOS_INLINE_FUNCTION
    double
    getImpurityRadius() const noexcept
    {
        return _impurityRadius;
    }

    virtual void
    setImpurityRadius(double impurityRadius) noexcept
    {
        _impurityRadius = impurityRadius;
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t
    getGridSize() const noexcept
    {
        return _gridSize;
    }

    virtual void
    setGridSize(std::size_t gridSize) = 0;

    virtual void
    setTemperatures(const std::vector<double>& gridTemperatures) = 0;

    virtual void
    syncClusterDataOnHost() = 0;

    virtual ClusterCommon<plsm::OnHost>
    getClusterCommon(std::size_t clusterId) const = 0;

    virtual ClusterCommon<plsm::OnHost>
    getSingleVacancy() = 0;

    virtual void
    computeAllFluxes(ConcentrationsView concentrations, FluxesView fluxes,
        std::size_t gridIndex) = 0;

    virtual void
    computeAllPartials(ConcentrationsView concentrations,
        Kokkos::View<double*> values, std::size_t gridIndex) = 0;

    virtual double
    getLargestRate() = 0;

    /**
     * Get the diagonal fill for the Jacobian, corresponding to the reactions.
     * Also populates the inverse map.
     *
     * @param fillMap Connectivity map.
     * @return The total number of partials.
     */
    virtual std::size_t
    getDiagonalFill(SparseFillMap& fillMap) = 0;

protected:
    double _latticeParameter {};
    double _atomicVolume {};
    double _interstitialBias {};
    double _impurityRadius {};

    std::size_t _gridSize {};
    std::size_t _numDOFs {};
    std::size_t _numClusters {};
};
}
}
