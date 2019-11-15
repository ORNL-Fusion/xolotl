#pragma once

#include <Kokkos_Core.hpp>

#include <plsm/detail/KokkosExtension.h>
#include <plsm/Utility.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TImpl>
template <typename TDerived>
class ReactionNetwork<TImpl>::Reaction
{
    static constexpr auto invalid = plsm::invalid<std::size_t>;

public:
    using NetworkType = ReactionNetwork<TImpl>;
    using Cluster = typename NetworkType::Cluster;
    using ConcentrationsView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
    using FluxesView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;

    enum class Type
    {
        production,
        dissociation
    };

    Reaction() = default;

    Reaction(NetworkType& network, std::size_t reactionId, Type reactionType,
        std::size_t cluster0, std::size_t cluster1,
        std::size_t cluster2 = invalid, std::size_t cluster3 = invalid);

    Type
    getType() const noexcept
    {
        return _type;
    }

    void
    updateRates();

    void
    productionFlux(ConcentrationsView concentrations, FluxesView fluxes)
    {
        using AmountType = typename TReactionNetwork::AmountType;
        constexpr auto numSpeciesNoI = TReactionNetwork::getNumberOfSpeciesNoI();
        constexpr auto speciesRangeNoI = TReactionNetwork::getSpeciesRangeNoI();

        // Compute the total number of elements in each cluster
        auto cl1 = _network.getCluster(_reactants[0]);
        const auto& cl1Reg = cl1.getRegion();
        AmountType nCl1 = 1;
        for (auto i : speciesRangeNoI) {
            nCl1 *= (cl1Reg[i].end() - 1 - cl1Reg[i].begin());
        }
        auto cl2 = _network.getCluster(_reactants[1]);
        const auto& cl2Reg = cl2.getRegion();
        AmountType nCl2 = 1;
        for (auto i : speciesRangeNoI) {
            nCl2 *= (cl2Reg[i].end() - 1 - cl2Reg[i].begin());
        }

        // Compute the flux for the 0th order moments
        double f = _coefs(0, 0, 0, 0) * concentrations[_reactants[0]] * concentrations[_reactants[1]];
        for (auto i : speciesRangeNoI) {
            f += _coefs(i() + 1, 0, 0, 0) * concentrations[_reactantMomentIds[0][i]] * concentrations[_reactants[1]];
        }
        for (auto j : speciesRangeNoI) {
            f += _coefs(0, j() + 1, 0, 0) * concentrations[_reactants[0]] * concentrations[_reactantMomentIds[1][j]];
        }
        for (auto i : speciesRangeNoI) {
            for (auto j : speciesRangeNoI) {
                f += _coefs(i() + 1, j() + 1, 0, 0) * concentrations[_reactantMomentIds[0][i]] * concentrations[_reactantMomentIds[1][j]];
            }
        }
        f *= _rate;

        fluxes[_reactants[0]] -= f / (double) nCl1;
        fluxes[_reactants[1]] -= f / (double) nCl2;
        for (auto prodId : _products) {
            if (prodId == invalid) {
                continue;
            }

            auto prod = _network.getCluster(prodId);
            const auto& prodReg = prod.getRegion();
            AmountType nProd = 1;
            for (auto i : speciesRangeNoI) {
                nProd *= (prodReg[i].end() - 1 - prodReg[i].begin());
            }

            fluxes[prodId] += f / (double) nProd;
        }

        // Take care of the first moments
        for (auto k : speciesRangeNoI) {
        	// First for the first reactant
            f = _coefs(0, 0, 1, k() + 1) * concentrations[_reactants[0]] * concentrations[_reactants[1]];
            for (auto i : speciesRangeNoI) {
                f += _coefs(i() + 1, 0, 1, k() + 1) * concentrations[_reactantMomentIds[0][i]] * concentrations[_reactants[1]];
            }
            for (auto j : speciesRangeNoI) {
                f += _coefs(0, j() + 1, 1, k() + 1) * concentrations[_reactants[0]] * concentrations[_reactantMomentIds[1][j]];
            }
            for (auto i : speciesRangeNoI) {
                for (auto j : speciesRangeNoI) {
                    f += _coefs(i() + 1, j() + 1, 1, k() + 1) * concentrations[_reactantMomentIds[0][i]] * concentrations[_reactantMomentIds[1][j]];
                }
            }
            f *= _rate;
            fluxes[_reactantMomentIds[0][k]] -= f / (double) nCl1;

        	// For the second reactant
            f = _coefs(0, 0, 2, k() + 1) * concentrations[_reactants[0]] * concentrations[_reactants[1]];
            for (auto i : speciesRangeNoI) {
                f += _coefs(i() + 1, 0, 2, k() + 1) * concentrations[_reactantMomentIds[0][i]] * concentrations[_reactants[1]];
            }
            for (auto j : speciesRangeNoI) {
                f += _coefs(0, j() + 1, 2, k() + 1) * concentrations[_reactants[0]] * concentrations[_reactantMomentIds[1][j]];
            }
            for (auto i : speciesRangeNoI) {
                for (auto j : speciesRangeNoI) {
                    f += _coefs(i() + 1, j() + 1, 2, k() + 1) * concentrations[_reactantMomentIds[0][i]] * concentrations[_reactantMomentIds[1][j]];
                }
            }
            f *= _rate;
            fluxes[_reactantMomentIds[1][k]] -= f / (double) nCl2;

        	// For the products
            for (auto prodId : _products) {
                if (prodId == invalid) {
                    continue;
                }

                auto prod = _network.getCluster(prodId);
                const auto& prodReg = prod.getRegion();
                AmountType nProd = 1;
                for (auto i : speciesRangeNoI) {
                    nProd *= (prodReg[i].end() - 1 - prodReg[i].begin());
                }

                // TODO generalize for 0, 1, or 2 products because here we are taking the same coefs for all of them
                f = _coefs(0, 0, 0, k() + 1) * concentrations[_reactants[0]] * concentrations[_reactants[1]];
                for (auto i : speciesRangeNoI) {
                    f += _coefs(i() + 1, 0, 0, k() + 1) * concentrations[_reactantMomentIds[0][i]] * concentrations[_reactants[1]];
                }
                for (auto j : speciesRangeNoI) {
                    f += _coefs(0, j() + 1, 0, k() + 1) * concentrations[_reactants[0]] * concentrations[_reactantMomentIds[1][j]];
                }
                for (auto i : speciesRangeNoI) {
                    for (auto j : speciesRangeNoI) {
                        f += _coefs(i() + 1, j() + 1, 0, k() + 1) * concentrations[_reactantMomentIds[0][i]] * concentrations[_reactantMomentIds[1][j]];
                    }
                }
                f *= _rate;
                // TODO use the correct first index
                fluxes[_productMomentIds[0][k]] -= f / (double) nProd;
            }
        }
    }

    void
    dissociationFlux(ConcentrationsView concentrations, FluxesView fluxes)
    {
        using AmountType = typename TReactionNetwork::AmountType;
        constexpr auto numSpeciesNoI = TReactionNetwork::getNumberOfSpeciesNoI();
        constexpr auto speciesRangeNoI = TReactionNetwork::getSpeciesRangeNoI();

        // Compute the total number of elements in each cluster
        auto cl = _network.getCluster(_reactants[0]);
        const auto& clReg = cl.getRegion();
        AmountType nCl = 1;
        for (auto i : speciesRangeNoI) {
            nCl *= (clReg[i].end() - 1 - clReg[i].begin());
        }
        auto prod1 = _network.getCluster(_products[0]);
        const auto& prod1Reg = prod1.getRegion();
        AmountType nProd1 = 1;
        for (auto i : speciesRangeNoI) {
            nProd1 *= (prod1Reg[i].end() - 1 - prod1Reg[i].begin());
        }
        auto prod2 = _network.getCluster(_products[1]);
        const auto& prod2Reg = prod2.getRegion();
        AmountType nProd2 = 1;
        for (auto i : speciesRangeNoI) {
            nProd2 *= (prod2Reg[i].end() - 1 - prod2Reg[i].begin());
        }

        // Compute the flux for the 0th order moments
        double f = _coefs(0, 0, 0, 0) * concentrations[_reactants[0]];
        for (auto i : speciesRangeNoI) {
            f += _coefs(i() + 1, 0, 0, 0) * concentrations[_reactantMomentIds[0][i]];
        }
        f *= _rate;
        fluxes[_reactants[0]] -= f / (double) nCl;
        fluxes[_products[0]] += f / (double) nProd1;
        fluxes[_products[1]] += f / (double) nProd2;

        // Take care of the first moments
        for (auto k : speciesRangeNoI) {
        	// First for the reactant
            f = _coefs(0, 0, 0, k() + 1) * concentrations[_reactants[0]];
            for (auto i : speciesRangeNoI) {
            	f += _coefs(i() + 1, 0, 0, k() + 1) * concentrations[_reactantMomentIds[0][i]];
            }
            // TODO compute the prefactor related to the dispersion, it can be moved to the coefs maybe
            double prefactor = 1.0;
            f *= _rate * prefactor;
            fluxes[_reactantMomentIds[0][k]] -= f / (double) nCl;

            // Now the first product
            f = _coefs(0, 0, 1, k() + 1) * concentrations[_reactants[0]];
            for (auto i : speciesRangeNoI) {
            	f += _coefs(i() + 1, 0, 1, k() + 1) * concentrations[_reactantMomentIds[0][i]];
            }
            // TODO compute the prefactor related to the dispersion, it can be moved to the coefs maybe
            prefactor = 1.0;
            f *= _rate * prefactor;
            fluxes[_productMomentIds[0][k]] += f / (double) nProd1;

            // Finally the second product
            f = _coefs(0, 0, 2, k() + 1) * concentrations[_reactants[0]];
            for (auto i : speciesRangeNoI) {
            	f += _coefs(i() + 1, 0, 2, k() + 1) * concentrations[_reactantMomentIds[0][i]];
            }
            // TODO compute the prefactor related to the dispersion, it can be moved to the coefs maybe
            prefactor = 1.0;
            f *= _rate * prefactor;
            fluxes[_productMomentIds[1][k]] += f / (double) nProd2;
        }
    }

    void
    contributeFlux(ConcentrationsView concentrations, FluxesView fluxes)
    {
        ((*this).*(_fluxFn))(concentrations, fluxes);
    }

private:
    TDerived*
    asDerived()
    {
        return static_cast<TDerived*>(this);
    }

    typename NetworkType::AmountType
    computeOverlap(Cluster singleCl, Cluster pairCl1, Cluster pairCl2);

    void
    computeProductionCoefficients();

    void
    computeDissociationCoefficients();

    void
    copyMomentIds(std::size_t clusterId,
        Kokkos::Array<std::size_t, 4>& momentIds)
    {
        if (clusterId == invalid) {
            momentIds = {invalid, invalid, invalid, invalid};
            return;
        }

        const auto& mIds = _network->getMomentIds(clusterId);
        for (std::size_t i = 0; i < 4; ++i) {
            momentIds[i] = mIds[i];
        }
    }

    double
    computeProductionRate(std::size_t gridIndex);

    double
    computeDissociationRate(std::size_t gridIndex);

    void
    computeProductionRates()
    {
        for (std::size_t i = 0; i < _rate.extent(0); ++i) {
            _rate(i) = asDerived()->computeProductionRate(i);
        }
    }

    void
    computeDissociationRates()
    {
        for (std::size_t i = 0; i < _rate.extent(0); ++i) {
            _rate(i) = asDerived()->computeDissociationRate(i);
        }
    }

private:
    NetworkType* _network {nullptr};

    Type _type {};

    using FluxFn = void (Reaction::*)(ConcentrationsView, FluxesView);
    FluxFn _fluxFn {nullptr};

    //Cluster indices for LHS and RHS of the reaction
    //Dissociation reactions always have 1 input and 2 outputs
    //Production reactions always have 2 inputs, but may have 0, 1, or 2 outputs
    Kokkos::Array<std::size_t, 2> _reactants {invalid, invalid};
    Kokkos::Array<std::size_t, 2> _products {invalid, invalid};

    Kokkos::Array<Kokkos::Array<std::size_t, 4>, 2> _reactantMomentIds;
    Kokkos::Array<Kokkos::Array<std::size_t, 4>, 2> _productMomentIds;

    //! Reaction rate (k)
    Kokkos::View<double*> _rate;

    //! Flux coefficients
    Kokkos::View<double****> _coefs;
};
}
}

#include <experimental/Reaction.inl>
