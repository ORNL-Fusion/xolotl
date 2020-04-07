#pragma once

#include <experimental/ReactionNetwork.h>
#include <experimental/FeReaction.h>
#include <experimental/FeTraits.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TSpeciesEnum>
class FeReactionGenerator;
}

template <typename TSpeciesEnum>
class FeReactionNetwork :
    public ReactionNetwork<FeReactionNetwork<TSpeciesEnum>>
{
    friend class ReactionNetwork<FeReactionNetwork<TSpeciesEnum>>;
    friend class detail::ReactionNetworkWorker<FeReactionNetwork<TSpeciesEnum>>;

public:
    using Superclass = ReactionNetwork<FeReactionNetwork<TSpeciesEnum>>;
    using Subpaving = typename Superclass::Subpaving;
    using Composition = typename Superclass::Composition;
    using Species = typename Superclass::Species;
    using AmountType = typename Superclass::AmountType;
    using IndexType = typename Superclass::IndexType;
    using ConcentrationsView = typename Superclass::ConcentrationsView;
    using FluxesView = typename Superclass::FluxesView;
    using ClusterDataRef = typename Superclass::ClusterDataRef;

    using Superclass::Superclass;

    void
    addReactionFluxes(ConcentrationsView concentrations, FluxesView fluxes,
        IndexType gridIndex)
    {
        auto nValues = this->_sinkingIds.extent(0);
        auto sinkingIds = this->_sinkingIds;
        auto sinkStrengths = this->_sinkStrengths;
        auto diffusionCoefficients = this->_clusterData.diffusionCoefficient;
        Kokkos::parallel_for(nValues, KOKKOS_LAMBDA (const IndexType i) {
            Kokkos::atomic_sub(&fluxes(sinkingIds(i)), concentrations(sinkingIds(i)) * sinkStrengths(i)
                * diffusionCoefficients(sinkingIds(i), gridIndex));
        });

        return;
    }

    void
    addReactionPartials(ConcentrationsView concentrations,
        Kokkos::View<double*> values, IndexType gridIndex)
    {
        auto nValues = this->_sinkingIds.extent(0);
        auto sinkingIds = this->_sinkingIds;
        auto sinkStrengths = this->_sinkStrengths;
        auto diffusionCoefficients = this->_clusterData.diffusionCoefficient;
        auto connectivity = this->_reactionData.connectivity;
        Kokkos::parallel_for(nValues, KOKKOS_LAMBDA (const IndexType i) {
            Kokkos::atomic_sub(&values(connectivity(sinkingIds(i), sinkingIds(i))), sinkStrengths(i)
                * diffusionCoefficients(sinkingIds(i), gridIndex));
        });
        return;
    }

private:
    double
    checkLatticeParameter(double latticeParameter)
    {
        if (latticeParameter <= 0.0) {
            return tungstenLatticeConstant;
        }
        return latticeParameter;
    }

    double
    checkImpurityRadius(double impurityRadius)
    {
        if (impurityRadius <= 0.0) {
            return heliumRadius;
        }
        return impurityRadius;
    }

    void
    initializeModifiedReactions()
    {
		double r0 = this->_latticeParameter * 0.75 * sqrt(3.0);
		double rho = 0.0003;
        // Initialize the views
        IndexType nSink = 5;
        _sinkingIds = Kokkos::View<IndexType*> ("sinkingIds", nSink);
        _sinkStrengths = Kokkos::View<double*> ("sinkStrengths", nSink);
        // Create mirror views
        auto hSinkingIds = create_mirror_view(_sinkingIds);
        auto hSinkStrengths = create_mirror_view(_sinkStrengths);

        // Look for the sinking clusters
        Composition comp = Composition::zero();
        comp[Species::I] = 1;
        auto cluster = this->findCluster(comp, plsm::onHost);
        hSinkingIds(0) = cluster.getId();
        hSinkStrengths(0) = 1.05 * -4.0 * xolotlCore::pi * rho
            / log(xolotlCore::pi * rho * (cluster.getReactionRadius() + r0)
            * (cluster.getReactionRadius() + r0));

        comp[Species::I] = 0;
        comp[Species::V] = 1;
        cluster = this->findCluster(comp, plsm::onHost);
        hSinkingIds(1) = cluster.getId();
        hSinkStrengths(1) = -4.0 * xolotlCore::pi * rho
                / log(xolotlCore::pi * rho * (cluster.getReactionRadius() + r0)
                * (cluster.getReactionRadius() + r0));

        comp[Species::V] = 2;
        cluster = this->findCluster(comp, plsm::onHost);
        hSinkingIds(2) = cluster.getId();
        hSinkStrengths(2) = -4.0 * xolotlCore::pi * rho
                / log(xolotlCore::pi * rho * (cluster.getReactionRadius() + r0)
                * (cluster.getReactionRadius() + r0));

        comp[Species::V] = 3;
        cluster = this->findCluster(comp, plsm::onHost);
        hSinkingIds(3) = cluster.getId();
        hSinkStrengths(3) = -4.0 * xolotlCore::pi * rho
                / log(xolotlCore::pi * rho * (cluster.getReactionRadius() + r0)
                * (cluster.getReactionRadius() + r0));

        comp[Species::V] = 4;
        cluster = this->findCluster(comp, plsm::onHost);
        hSinkingIds(4) = cluster.getId();
        hSinkStrengths(4) = -4.0 * xolotlCore::pi * rho
                / log(xolotlCore::pi * rho * (cluster.getReactionRadius() + r0)
                * (cluster.getReactionRadius() + r0));

        deep_copy(_sinkingIds, hSinkingIds);
        deep_copy(_sinkStrengths, hSinkStrengths);

        return;
    }

    detail::FeReactionGenerator<Species>
    getReactionGenerator() const noexcept
    {
        return detail::FeReactionGenerator<Species>{*this};
    }

public:

    Kokkos::View<IndexType*> _sinkingIds;
    Kokkos::View<double*> _sinkStrengths;
};

namespace detail
{
template <typename TSpeciesEnum>
class FeReactionGenerator : public
    ReactionGenerator<FeReactionNetwork<TSpeciesEnum>,
        FeReactionGenerator<TSpeciesEnum>>
{
public:
    using Network = FeReactionNetwork<TSpeciesEnum>;
    using Subpaving = typename Network::Subpaving;
    using IndexType = typename Network::IndexType;

    using Superclass = ReactionGenerator<FeReactionNetwork<TSpeciesEnum>,
        FeReactionGenerator<TSpeciesEnum>>;

    using Superclass::Superclass;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    operator()(IndexType i, IndexType j, TTag tag) const;
};
}
}
}

#include <experimental/FeClusterGenerator.h>
#include <experimental/FeReactionNetwork.inl>
