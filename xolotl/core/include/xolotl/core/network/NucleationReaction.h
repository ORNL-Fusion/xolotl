#pragma once

#include <xolotl/core/network/Reaction.h>
#include <xolotl/core/network/SpeciesEnumSequence.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
class NucleationReaction : public Reaction<TNetwork, TDerived>
{
    friend class Reaction<TNetwork, TDerived>;

public:
    using NetworkType = TNetwork;
    using Superclass = Reaction<TNetwork, TDerived>;
    using ClusterDataRef = typename Superclass::ClusterDataRef;
    using IndexType = typename Superclass::IndexType;
    using Connectivity = typename Superclass::Connectivity;
    using ConcentrationsView = typename Superclass::ConcentrationsView;
    using FluxesView = typename Superclass::FluxesView;
    using AmountType = typename Superclass::AmountType;
    using ReactionDataRef = typename Superclass::ReactionDataRef;

    NucleationReaction() = default;

    KOKKOS_INLINE_FUNCTION
    NucleationReaction(ReactionDataRef reactionData,
            ClusterDataRef clusterData, IndexType reactionId,
            IndexType cluster0,
            IndexType cluster1)
        :
        Superclass(reactionData, clusterData, reactionId),
        _reactant(cluster0),
        _product(cluster1)
    {
        this->initialize();
    }

    KOKKOS_INLINE_FUNCTION
    NucleationReaction(ReactionDataRef reactionData,
            ClusterDataRef clusterData, IndexType reactionId,
            const detail::ClusterSet& clusterSet)
        :
        NucleationReaction(reactionData, clusterData, reactionId, clusterSet.cluster0, clusterSet.cluster1)
    {
    }

    static
    detail::CoefficientsView
    allocateCoefficientsView(IndexType)
    {
        return detail::CoefficientsView();
    }

private:
    KOKKOS_INLINE_FUNCTION
    void
    computeCoefficients()
    {
        // No coefs
    }

    KOKKOS_INLINE_FUNCTION
    double
    computeRate(IndexType gridIndex)
    {
        // We say there are 25 bubbles created per fission fragments and there are 2 fission fragments per fission
        double rate = 50.0 * this->_clusterData.fissionRate(0);

        return rate;
    }

    KOKKOS_INLINE_FUNCTION
    void
    computeConnectivity(const Connectivity& connectivity)
    {
        // The reactant connects with the reactant
        this->addConnectivity(_reactant, _reactant, connectivity);
        // The product connects with the reactant
        this->addConnectivity(_product, _reactant, connectivity);
    }

    KOKKOS_INLINE_FUNCTION
    void
    computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
        IndexType gridIndex)
    {
        // Get the single concentration to know in which regime we are
        double singleConc = concentrations(_reactant);

    	// Update the concentrations
    	if (singleConc > 2.0 * this->_rate(gridIndex)) {
    		Kokkos::atomic_sub(&fluxes(_reactant), 2.0 * this->_rate(gridIndex));
    		Kokkos::atomic_add(&fluxes(_product), this->_rate(gridIndex));
    	} else {
    		Kokkos::atomic_sub(&fluxes(_reactant), singleConc);
    		Kokkos::atomic_add(&fluxes(_product), singleConc / 2.0);
    	}
    }

    KOKKOS_INLINE_FUNCTION
    void
    computePartialDerivatives(ConcentrationsView concentrations,
        Kokkos::View<double*> values, Connectivity connectivity,
        IndexType gridIndex)
    {

    	// Get the single concentration to know in which regime we are
    	double singleConc = concentrations(_reactant);

    	// Update the partials
    	if (singleConc > 2.0 * this->_rate(gridIndex)) {
            // Nothing
    	} else {
    		Kokkos::atomic_sub(&values(connectivity(_reactant, _reactant)), 1.0);
    		Kokkos::atomic_add(&values(connectivity(_product, _reactant)), 0.5);
    	}
    }

    KOKKOS_INLINE_FUNCTION
    double
    computeLeftSideRate(ConcentrationsView concentrations, IndexType clusterId,
        IndexType gridIndex)
    {
        return 0.0;
    }

protected:
    IndexType _reactant;
    IndexType _product;
    static constexpr auto invalidIndex = Superclass::invalidIndex;
};
}
}
}

#include <xolotl/core/network/detail/NucleationReactionGenerator.h>
