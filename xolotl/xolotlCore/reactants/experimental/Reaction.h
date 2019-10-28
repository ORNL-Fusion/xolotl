#pragma once


namespace xolotlCore
{
namespace experimental
{
template <typename TDerived>
class Reaction
{
    static constexpr auto invalid = plsm::invalid<std::size_t>;

public:
    using ConcentrationsView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
    using FluxesView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;

    enum class Type
    {
        production,
        dissociation
    };

    Reaction() = default;

    template <typename TReactionNetwork>
    Reaction(TReactionNetwork& network, Type reactionType,
        std::size_t cluster0, std::size_t cluster1,
        std::size_t cluster2 = invalid, std::size_t cluster3 = invalid);

    template <typename TReactionNetwork>
    void
    computeProductionCoefficients(TReactionNetwork& network);

    template <typename TReactionNetwork>
    void
    computeDissociationCoefficients(TReactionNetwork& network);

    Type
    getType() const noexcept
    {
        return _type;
    }

    void
    productionFlux(ConcentrationsView concentrations, FluxesView fluxes)
    {
        //TODO: Compute flux
        double f = 1.0;
        fluxes[_reactants[0]] -= f;
        fluxes[_reactants[1]] -= f;
        for (auto prodId : _products) {
            if (prodId == invalid) {
                continue;
            }
            fluxes[prodId] += f;
        }
    }

    void
    dissociationFlux(ConcentrationsView concentrations, FluxesView fluxes)
    {
        double f = 1.0;
        fluxes[_reactants[0]] -= f;
        fluxes[_products[0]] += f;
        fluxes[_products[1]] += f;
    }

    void
    contributeFlux(ConcentrationsView concentrations, FluxesView fluxes)
    {
        ((*this).*(_fluxFn))(concentrations, fluxes);
    }

private:
    Type _type {};

    using FluxFn = void (Reaction::*)(ConcentrationsView, FluxesView);
    FluxFn _fluxFn {nullptr};

    //Cluster indices for LHS and RHS
    //Dissociation reactions always have 1 input and 2 outputs
    //Production reactions always have 2 inputs, but may have 0, 1, or 2 outputs
    Kokkos::Array<std::size_t, 2> _reactants {invalid, invalid};
    Kokkos::Array<std::size_t, 2> _products {invalid, invalid};

    //! Reaction rate, k
    double _rate {};

    //! Flux coefficients
    Kokkos::View<double****> _coefs;
};


class PSIReaction : public Reaction<PSIReaction>
{
public:
    using Reaction<PSIReaction>::Reaction;

    template <typename TReactionNetwork>
    double
    computeRate(TReactionNetwork& network)
    {
        return 1.0;
    }
};
}
}

#include <experimental/Reaction.inl>
