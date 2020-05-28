#pragma once

#include <xolotl/core/reactants/ReactionNetwork.h>
#include <xolotl/core/reactants/NEReaction.h>
#include <xolotl/core/reactants/NETraits.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
class NEReactionGenerator;

class NEClusterUpdater;
}

class NEReactionNetwork : public ReactionNetwork<NEReactionNetwork>
{
    friend class ReactionNetwork<NEReactionNetwork>;
    friend class detail::ReactionNetworkWorker<NEReactionNetwork>;

public:
    using Superclass = ReactionNetwork<NEReactionNetwork>;
    using Subpaving = typename Superclass::Subpaving;
    using Composition = typename Superclass::Composition;
    using Species = typename Superclass::Species;
    using IndexType = typename Superclass::IndexType;
    using ConcentrationsView = typename Superclass::ConcentrationsView;
    using FluxesView = typename Superclass::FluxesView;

    using Superclass::Superclass;

private:
    double
    checkLatticeParameter(double latticeParameter)
    {
        if (latticeParameter <= 0.0) {
            return uraniumDioxydeLatticeConstant;
        }
        return latticeParameter;
    }

    double checkImpurityRadius(double impurityRadius)
    {
        if (impurityRadius <= 0.0) {
            return xenonRadius;
        }
        return impurityRadius;
    }

    void
    checkTiles(const IOptions& options)
    {
        return;
    }

    detail::NEReactionGenerator
    getReactionGenerator() const noexcept;
};

namespace detail
{
class NEReactionGenerator :
    public ReactionGenerator<NEReactionNetwork, NEReactionGenerator>
{
    friend class ReactionGeneratorBase<NEReactionNetwork, NEReactionGenerator>;

public:
    using NetworkType = NEReactionNetwork;
    using Subpaving = typename NetworkType::Subpaving;
    using Superclass =
        ReactionGenerator<NEReactionNetwork, NEReactionGenerator>;

    using Superclass::Superclass;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    operator()(IndexType i, IndexType j, TTag tag) const;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    addSinks(IndexType i, TTag tag) const;

private:
    ReactionCollection<NetworkType>
    getReactionCollection() const;
};

class NEClusterUpdater
{
public:
    using NetworkType = NEReactionNetwork;
    using ClusterData = typename NetworkType::ClusterData;
    using IndexType = typename NetworkType::IndexType;

    KOKKOS_INLINE_FUNCTION
    void
    updateDiffusionCoefficient(const ClusterData& data, IndexType clusterId,
        IndexType gridIndex) const
    {
        // If the diffusivity is given
        if (data.migrationEnergy(clusterId) > 0.0) {

            // Intrinsic diffusion
            double kernel = -3.04 / (kBoltzmann * data.temperature(gridIndex));
            double D3 = 7.6e8 * exp(kernel); // nm2/s

            // We need the fission rate now
            double fissionRate = data.fissionRate(0) * 1.0e27; // #/m3/s

            // Athermal diffusion
            double D1 = (8e-40 * fissionRate) * 1.0e18; // nm2/s

            // Radiation-enhanced diffusion
            kernel = -1.2 / (kBoltzmann * data.temperature(gridIndex));
            double D2 = (5.6e-25 * sqrt(fissionRate) * exp(kernel)) * 1.0e18; // nm2/s

            data.diffusionCoefficient(clusterId, gridIndex) = D1 + D2 + D3;

            return;
        }

        data.diffusionCoefficient(clusterId, gridIndex) = data.diffusionFactor(clusterId);
    }
};
}
}
}

#include <xolotl/core/reactants/NEClusterGenerator.h>
#include <xolotl/core/reactants/NEReactionNetwork.inl>
