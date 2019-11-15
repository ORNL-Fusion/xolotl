#pragma once

#include <algorithm>
#include <array>

#include <plsm/EnumIndexed.h>

#include <Constants.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TImpl>
template <typename TDerived>
ReactionNetwork<TImpl>::Reaction<TDerived>::Reaction(NetworkType& network,
    std::size_t reactionId, Type reactionType, std::size_t cluster0,
    std::size_t cluster1, std::size_t cluster2, std::size_t cluster3)
    :
    _network(&network),
    _type(reactionType),
    _fluxFn(
        _type == Type::production ? &Reaction::productionFlux :
            &Reaction::dissociationFlux),
    _reactants(
        _type == Type::production ? Kokkos::Array<std::size_t, 2>( {cluster0,
            cluster1}) : Kokkos::Array<std::size_t, 2>( {cluster0, invalid})),
    _products(_type == Type::production ? Kokkos::Array<std::size_t, 2>( {
        cluster2, cluster3}) : Kokkos::Array<std::size_t, 2>( {cluster1,
        cluster2})),
    _rate(_network->getReactionRates(reactionId))
{
    for (std::size_t i : {0, 1}) {
        copyMomentIds(_reactants[i], _reactantMomentIds[i]);
        copyMomentIds(_products[i], _productMomentIds[i]);
    }

    //TODO: Replace '5' with (numSpeciesNoI + 1)
    if (_type == Type::production) {
        int nProd = 0;
        for (auto prodId : _products) {
            if (prodId != invalid) {
                ++nProd;
            }
        }
        auto nCl = 2 + nProd;
        _coefs = Kokkos::View<double****>("Flux Coefficients", 5, 5, nCl, 5);
        computeProductionCoefficients();
    }
    else {
        _coefs = Kokkos::View<double****>("Flux Coefficients", 5, 1, 3, 5);
        computeDissociationCoefficients();
    }

    updateRates();
}

template <typename TImpl>
template <typename TDerived>
inline void
ReactionNetwork<TImpl>::Reaction<TDerived>::updateRates()
{
    if (_type == Type::production) {
        computeProductionRates();
    }
    else {
        computeDissociationRates();
    }
}

template <typename TImpl>
template <typename TDerived>
inline typename ReactionNetwork<TImpl>::AmountType
ReactionNetwork<TImpl>::Reaction<TDerived>::computeOverlap(
    Cluster singleCl, Cluster pairCl1, Cluster pairCl2)
{
    using AmountType = typename NetworkType::AmountType;
    constexpr auto numSpeciesNoI = NetworkType::getNumberOfSpeciesNoI();
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

    const auto& sReg = singleCl.getRegion();
    const auto& pReg1 = pairCl1.getRegion();
    const auto& pReg2 = pairCl2.getRegion();

    AmountType nOverlap = 1;
    for (auto i : speciesRangeNoI) {
        // The width is the subset of the tiles for which the
        // reaction is possible
        // For instance if we have X_1 + X_[3,5) ⇄ X_[5,7)
        // It is only possible for 4 within X_[3,5) and 5 within X_[5,7)
        // so the width is 1
        // More complicated with X_[3,5) + X_[5,7) ⇄ X_[9,11)
        // 3+6, 4+5, 4+6, width is 3

        AmountType width{};

        //TODO: Would be nice to loop on the cluster with the smaller tile
        for (auto j : makeIntervalRange(pReg1[i])) {
            width += std::min(sReg[i].end() - 1, pReg2[i].end() - 1 + j)
                - std::max(sReg[i].begin(), pReg2[i].begin() + j) + 1;
        }

        nOverlap *= width;
    }

    assert(nOverlap > 0);

    return nOverlap;
}

template <typename TImpl>
template <typename TDerived>
inline void
ReactionNetwork<TImpl>::Reaction<TDerived>::computeProductionCoefficients()
{
    using Species = typename NetworkType::Species;
    using SpeciesSeq = typename NetworkType::SpeciesSequence;
    using AmountType = typename NetworkType::AmountType;

    // Find the overlap for this reaction
    constexpr auto numSpeciesNoI = NetworkType::getNumberOfSpeciesNoI();
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

    auto cl1 = _network->getCluster(_reactants[0]);
    auto cl2 = _network->getCluster(_reactants[1]);
    auto prod = _network->getCluster(_products[0]);

    auto nOverlap = static_cast<double>(computeOverlap(prod, cl1, cl2));

    const auto& cl1Reg = cl1.getRegion();
    const auto& cl2Reg = cl2.getRegion();
    const auto& prodReg = prod.getRegion();

    _coefs(0, 0, 0, 0) = nOverlap;
    for (auto i : speciesRangeNoI) {
        // First order sum
        for (auto l : makeIntervalRange(cl2Reg[i])) {
            _coefs(i() + 1, 0, 0, 0) += firstOrderSum(
                std::max(prodReg[i].begin() - l, cl1Reg[i].begin()),
                std::min(prodReg[i].end() - 1 - l, cl1Reg[i].end() - 1),
                static_cast<double>(cl1Reg[i].end() - 1 + cl1Reg[i].begin())
                    / 2.0);
        }

        for (auto l : makeIntervalRange(cl1Reg[i])) {
            _coefs(0, i() + 1, 0, 0) += firstOrderSum(
                std::max(prodReg[i].begin() - l, cl2Reg[i].begin()),
                std::min(prodReg[i].end() - 1 - l, cl2Reg[i].end() - 1),
                static_cast<double>(cl2Reg[i].end() - 1 + cl2Reg[i].begin())
                    / 2.0);
        }

        for (auto l : makeIntervalRange(cl1Reg[i])) {
            _coefs(0, 0, 0, i() + 1) += firstOrderSum(
                std::max(prodReg[i].begin(), cl2Reg[i].begin() + l),
                std::min(prodReg[i].end() - 1, cl2Reg[i].end() - 1 + l),
                static_cast<double>(prodReg[i].end() - 1 + prodReg[i].begin())
                    / 2.0);
        }

        // Product partial derivatives
        for (auto k : speciesRangeNoI) {
            // Second order sum
            if (k == i) {
                for (auto l : makeIntervalRange(cl2Reg[i])) {
                    _coefs(i() + 1, 0, 0, k() + 1) += secondOrderOffsetSum(
                        std::max(prodReg[i].begin() - l, cl1Reg[i].begin()),
                        std::min(prodReg[i].end() - 1 - l, cl1Reg[i].end() - 1),
                        (double) (cl1Reg[i].end() - 1 + cl1Reg[i].begin())
                            / 2.0,
                        (double) (prodReg[i].end() - 1 + prodReg[i].begin())
                            / 2.0, l);
                }
                for (auto l : makeIntervalRange(cl1Reg[i])) {
                    _coefs(0, i() + 1, 0, k() + 1) += secondOrderOffsetSum(
                        std::max(prodReg[i].begin() - l, cl2Reg[i].begin()),
                        std::min(prodReg[i].end() - 1 - l, cl2Reg[i].end() - 1),
                        (double) (cl2Reg[i].end() - 1 + cl2Reg[i].begin())
                            / 2.0,
                        (double) (prodReg[i].end() - 1 + prodReg[i].begin())
                            / 2.0, l);
                }
            }
            else {
                _coefs(i() + 1, 0, 0, k() + 1) += _coefs(i() + 1, 0, 0, 0)
                    * _coefs(0, 0, 0, k() + 1) / nOverlap;

                _coefs(0, i() + 1, 0, k() + 1) += _coefs(0, i() + 1, 0, 0)
                    * _coefs(0, 0, 0, k() + 1) / nOverlap;
            }
        }
    }

    for (auto i : speciesRangeNoI) {
        // First reactant partial derivatives
        for (auto k : speciesRangeNoI) {
            _coefs(0, 0, 1, k() + 1) += _coefs(k() + 1, 0, 0, 0);

            if (k == i) {
                for (auto l : makeIntervalRange(cl2Reg[i])) {
                    _coefs(i() + 1, 0, 1, k()) += secondOrderSum(
                        std::max(prodReg[i].begin() - l, cl1Reg[i].begin()),
                        std::min(prodReg[i].end() - 1 - l, cl1Reg[i].end() - 1),
                        (double) (cl1Reg[i].end() - 1 + cl1Reg[i].begin())
                            / 2.0);
                }
            }
            else {
                _coefs(i() + 1, 0, 1, k() + 1) += _coefs(i() + 1, 0, 0, 0)
                    * _coefs(k() + 1, 0, 0, 0) / nOverlap;
            }

            _coefs(0, i() + 1, 1, k() + 1) += _coefs(k() + 1, i(), 0, 0);
        }

        // Second reactant partial derivatives
        for (auto k : speciesRangeNoI) {
            _coefs(0, 0, 2, k() + 1) += _coefs(0, k() + 1, 0, 0);

            if (k == i) {
                for (auto l : makeIntervalRange(cl1Reg[i])) {
                    _coefs(0, i() + 1, 2, k() + 1) += secondOrderSum(
                        std::max(prodReg[i].begin() - l, cl2Reg[i].begin()),
                        std::min(prodReg[i].end() - 1 - l, cl2Reg[i].end() - 1),
                        (double) (cl2Reg[i].end() - 1 + cl2Reg[i].begin())
                            / 2.0);
                }
            }
            else {
                _coefs(0, i() + 1, 2, k() + 1) += _coefs(0, i() + 1, 0, 0)
                    * _coefs(0, k() + 1, 0, 0) / nOverlap;
            }

            _coefs(i() + 1, 0, 2, k() + 1) += _coefs(i() + 1, k() + 1, 0, 0);
        }
    }

    // Now we loop over the 2 dimensions of the coefs to compute all
    // the possible sums over distances for the flux
    // The first index represents the first reactant,
    // the second one is the second reactant
    // 0 is the 0th order, 1, 2, 3, 4, are the distances in the
    // He, D, T, V directions (TODO make sure the order is correct)
    for (auto i : speciesRangeNoI) {
        for (auto j : speciesRangeNoI) {
            // Second order sum
            if (i == j) {
                for (auto l : makeIntervalRange(cl1Reg[j])) {
                    _coefs(i() + 1, j() + 1, 0, 0) += (l
                        - (double) (cl1Reg[j].end() - 1 + cl1Reg[j].begin())
                            / 2.0)
                        * firstOrderSum(
                            std::max(prodReg[j].begin() - l, cl2Reg[j].begin()),
                            std::min(prodReg[j].end() - 1 - l,
                                cl2Reg[j].end() - 1),
                            (double) (cl2Reg[j].end() - 1 + cl2Reg[j].begin())
                                / 2.0);
                }
            }
            else {
                _coefs(i() + 1, j() + 1, 0, 0) += _coefs(i() + 1, 0, 0, 0)
                    * _coefs(0, j() + 1, 0, 0) / nOverlap;
            }

            // Now we deal with the coefficients needed for the
            // partial derivatives
            // Let's start with the product
            for (auto k : speciesRangeNoI) {
                // Third order sum
                if (i == j && j == k) {
                    for (auto l : makeIntervalRange(cl1Reg[i])) {
                        _coefs(i() + 1, j() + 1, 0, k() + 1) += (l
                            - (double) (cl1Reg[i].end() - 1 + cl1Reg[i].begin())
                                / 2.0)
                            * secondOrderOffsetSum(
                                std::max(prodReg[i].begin() - l,
                                    cl2Reg[i].begin()),
                                std::min(prodReg[i].end() - 1 - l,
                                    cl2Reg[i].end() - 1),
                                (double) (cl2Reg[i].end() - 1
                                    + cl2Reg[i].begin()) / 2.0,
                                (double) (prodReg[i].end() - 1
                                    + prodReg[i].begin()) / 2.0, l);
                    }
                }
                else if (j == k) {
                    _coefs(i() + 1, j() + 1, 0, k() + 1) += _coefs(i() + 1, 0,
                        0, 0) * _coefs(0, j() + 1, 0, k() + 1) / nOverlap;
                }
                else if (i == k) {
                    _coefs(i() + 1, j() + 1, 0, k() + 1) += _coefs(0, j() + 1,
                        0, 0) * _coefs(i() + 1, 0, 0, k() + 1) / nOverlap;
                }
                else {
                    // TODO check this is the right formula, might be divided by nOverlap^2
                    _coefs(i() + 1, j() + 1, 0, k() + 1) += _coefs(i() + 1, 0,
                        0, 0) * _coefs(0, j() + 1, 0, 0)
                        * _coefs(0, 0, 0, k() + 1) / nOverlap;
                }
            }
            // Let's take care of the first reactant partial derivatives
            for (auto k : speciesRangeNoI) {
                // Third order sum
                if (i == j && j == k) {
                    for (auto l : makeIntervalRange(cl1Reg[i])) {
                        _coefs(i() + 1, j() + 1, 1, k() + 1) += (l
                            - (double) (cl1Reg[i].end() - 1 + cl1Reg[i].begin())
                                / 2.0)
                            * (l
                                - (double) (cl1Reg[i].end() - 1
                                    + cl1Reg[i].begin()) / 2.0)
                            * firstOrderSum(
                                std::max(prodReg[i].begin() - l,
                                    cl2Reg[i].begin()),
                                std::min(prodReg[i].end() - 1 - l,
                                    cl2Reg[i].end() - 1),
                                (double) (cl2Reg[i].end() - 1
                                    + cl2Reg[i].begin()) / 2.0);
                    }
                }
                else if (i == k) {
                    _coefs(i() + 1, j() + 1, 1, k() + 1) += _coefs(0, j() + 1,
                        0, 0) * _coefs(i() + 1, 0, 1, k() + 1) / nOverlap;
                }
                else if (j == k) {
                    _coefs(i() + 1, j() + 1, 1, k() + 1) += _coefs(i() + 1, 0,
                        0, 0) * _coefs(0, j() + 1, 1, k() + 1) / nOverlap;
                }
                else {
                    // TODO check this is the right formula, might be divided by nOverlap^2
                    _coefs(i() + 1, j() + 1, 1, k() + 1) += _coefs(i() + 1, 0,
                        0, 0) * _coefs(0, j() + 1, 0, 0)
                        * _coefs(k() + 1, 0, 0, 0) / nOverlap;
                }
            }

            // Let's take care of the second reactant partial derivatives
            for (auto k : speciesRangeNoI) {
                // Third order sum
                if (i == j && j == k) {
                    for (auto l : makeIntervalRange(cl2Reg[i])) {
                        _coefs(i() + 1, j() + 1, 2, k() + 1) += (l
                            - (double) (cl2Reg[i].end() - 1 + cl2Reg[i].begin())
                                / 2.0)
                            * (l
                                - (double) (cl2Reg[i].end() - 1
                                    + cl2Reg[i].begin()) / 2.0)
                            * firstOrderSum(
                                std::max(prodReg[i].begin() - l,
                                    cl1Reg[i].begin()),
                                std::min(prodReg[i].end() - 1 - l,
                                    cl1Reg[i].end() - 1),
                                (double) (cl1Reg[i].end() - 1
                                    + cl1Reg[i].begin()) / 2.0);
                    }
                }
                else if (i == k) {
                    _coefs(i() + 1, j() + 1, 2, k() + 1) += _coefs(0, j() + 1,
                        0, 0) * _coefs(i() + 1, 0, 2, k() + 1) / nOverlap;
                }
                else if (j == k) {
                    _coefs(i() + 1, j() + 1, 2, k() + 1) += _coefs(i() + 1, 0,
                        0, 0) * _coefs(0, j() + 1, 2, k() + 1) / nOverlap;
                }
                else {
                    // TODO check this is the right formula, might be divided by nOverlap^2
                    _coefs(i() + 1, j() + 1, 2, k() + 1) += _coefs(i() + 1, 0,
                        0, 0) * _coefs(0, j() + 1, 0, 0)
                        * _coefs(0, k() + 1, 0, 0) / nOverlap;
                }
            }
        }
    }
}

template <typename TImpl>
template <typename TDerived>
inline void
ReactionNetwork<TImpl>::Reaction<TDerived>::computeDissociationCoefficients()
{
    using Species = typename NetworkType::Species;
    using SpeciesSeq = typename NetworkType::SpeciesSequence;
    using AmountType = typename NetworkType::AmountType;

    constexpr auto numSpeciesNoI = NetworkType::getNumberOfSpeciesNoI();
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

    auto cl = _network->getCluster(_reactants[0]);
    auto prod1 = _network->getCluster(_products[0]);
    auto prod2 = _network->getCluster(_products[1]);

    auto nOverlap = static_cast<double>(computeOverlap(cl, prod1, prod2));

    auto clReg = cl.getRegion();
    auto prod1Reg = prod1.getRegion();
    auto prod2Reg = prod2.getRegion();

    // The first coefficient is simply the overlap because it is the sum over 1
    _coefs(0, 0, 0, 0) = nOverlap;
    for (auto i : speciesRangeNoI) {
        // First order sum
        for (auto l : makeIntervalRange(prod1Reg[i])) {
            _coefs(i() + 1, 0, 0, 0) += firstOrderSum(
                std::max(clReg[i - 1].begin(), prod2Reg[i - 1].begin() + l),
                std::min(clReg[i - 1].end() - 1, prod2Reg[i - 1].end() - 1 + l),
                (double) (clReg[i - 1].end() - 1 + clReg[i - 1].begin()) / 2.0);
        }
    }

    // Partial derivatives
    for (auto k : speciesRangeNoI) {
        // Reactant
        _coefs(0, 0, 0, k() + 1) += _coefs(k() + 1, 0, 0, 0);

        // First product
        for (auto l : makeIntervalRange(prod2Reg[k])) {
            _coefs(0, 0, 1, k() + 1) += firstOrderSum(
                std::max(clReg[k].begin() - l, prod1Reg[k].begin()),
                std::min(clReg[k].end() - 1 - l, prod1Reg[k].end() - 1),
                (double) (prod1Reg[k].end() - 1 + prod1Reg[k].begin()) / 2.0);
        }

        // Second product
        for (auto l : makeIntervalRange(prod1Reg[k])) {
            _coefs(0, 0, 2, k() + 1) += firstOrderSum(
                std::max(clReg[k].begin() - l, prod2Reg[k].begin()),
                std::min(clReg[k].end() - 1 - l, prod2Reg[k].end() - 1),
                (double) (prod2Reg[k].end() - 1 + prod2Reg[k].begin()) / 2.0);
        }
    }

    // Now we loop over the 1 dimension of the coefs to compute all the
    // possible sums over distances for the flux
    // The first index represents the reactant
    // 0 is the 0th order, 1, 2, 3, 4, are the distances in the
    // He, D, T, V directions
    // TODO: make sure the order is correct
    for (auto i : speciesRangeNoI) {
        // Now we deal with the coefficients needed for the partial derivatives
        // Starting with the reactant
        for (auto k : speciesRangeNoI) {
            // Second order sum
            if (k == i) {
                for (auto l : makeIntervalRange(prod1Reg[i])) {
                    _coefs(i() + 1, 0, 0, k() + 1) += secondOrderSum(
                        std::max(clReg[i].begin(), prod2Reg[i].begin() + l),
                        std::min(clReg[i].end() - 1, prod2Reg[i].end() - 1 + l),
                        (double) (clReg[i].end() - 1 + clReg[i].begin()) / 2.0);
                }
            }
            else {
                _coefs(i() + 1, 0, 0, k() + 1) += _coefs(i() + 1, 0, 0, 0)
                    * _coefs(k() + 1, 0, 0, 0) / nOverlap;
            }
        }

        // Partial derivatives for the first product
        for (auto k : speciesRangeNoI) {
            // Second order sum
            if (k == i) {
                for (auto l : makeIntervalRange(prod2Reg[i])) {
                    _coefs(i() + 1, 0, 1, k() + 1) += secondOrderOffsetSum(
                        std::max(clReg[i].begin(), prod1Reg[i].begin() + l),
                        std::min(clReg[i].end() - 1, prod1Reg[i].end() - 1 + l),
                        (double) (clReg[i].end() - 1 + clReg[i].begin()) / 2.0,
                        (double) (prod1Reg[i].end() - 1 + prod1Reg[i].begin())
                            / 2.0, -l);
                }
            }
            else {
                _coefs(i() + 1, 0, 1, k() + 1) += _coefs(i() + 1, 0, 0, 0)
                    * _coefs(0, 0, 1, k() + 1) / nOverlap;
            }
        }

        // Partial derivatives for the second product
        for (auto k : speciesRangeNoI) {
            // Second order sum
            if (k == i) {
                for (auto l : makeIntervalRange(prod1Reg[i])) {
                    _coefs(i() + 1, 0, 2, k() + 1) += secondOrderOffsetSum(
                        std::max(clReg[i].begin(), prod2Reg[i].begin() + l),
                        std::min(clReg[i].end() - 1, prod2Reg[i].end() - 1 + l),
                        (double) (clReg[i].end() - 1 + clReg[i].begin()) / 2.0,
                        (double) (prod2Reg[i].end() - 1 + prod2Reg[i].begin())
                            / 2.0, -l);
                }
            }
            else {
                _coefs(i() + 1, 0, 2, k() + 1) += _coefs(i() + 1, 0, 0, 0)
                    * _coefs(0, 0, 2, k() + 1) / nOverlap;
            }
        }
    }
}

template <typename TImpl>
template <typename TDerived>
inline double
ReactionNetwork<TImpl>::Reaction<TDerived>::computeProductionRate(
    std::size_t gridIndex)
{
    double r0 = _network->getReactionRadius(_reactants[0]);
    double r1 = _network->getReactionRadius(_reactants[1]);

    double dc0 = _network->getDiffusionCoefficient(_reactants[0], gridIndex);
    double dc1 = _network->getDiffusionCoefficient(_reactants[1], gridIndex);

    constexpr double pi = ::xolotlCore::pi;

    double kPlus = 4.0 * pi * (r0 + r1) * (dc0 + dc1);

    return kPlus;
}

template <typename TImpl>
template <typename TDerived>
inline double
ReactionNetwork<TImpl>::Reaction<TDerived>::computeDissociationRate(
    std::size_t gridIndex)
{
    double omega = _network->getAtomicVolume();
    double T = _network->getTemperature(gridIndex);

    double kPlus = asDerived()->computeProductionRate(gridIndex);
    double E_b = asDerived()->computeBindingEnergy();

    constexpr double k_B = ::xolotlCore::kBoltzmann;

    double kMinus = (1.0 / omega) * kPlus * std::exp(-E_b / (k_B * T));

    return kMinus;
}

template <typename TImpl>
template <typename TDerived>
inline void
ReactionNetwork<TImpl>::Reaction<TDerived>::productionFlux(
    ConcentrationsView concentrations, FluxesView fluxes, std::size_t gridIndex)
{
    using AmountType = typename NetworkType::AmountType;
    constexpr auto numSpeciesNoI = NetworkType::getNumberOfSpeciesNoI();
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

    // Compute the total number of elements in each cluster
    auto cl1 = _network->getCluster(_reactants[0]);
    const auto& cl1Reg = cl1.getRegion();
    AmountType nCl1 = 1;
    for (auto i : speciesRangeNoI) {
        nCl1 *= (cl1Reg[i].end() - 1 - cl1Reg[i].begin());
    }
    auto cl2 = _network->getCluster(_reactants[1]);
    const auto& cl2Reg = cl2.getRegion();
    AmountType nCl2 = 1;
    for (auto i : speciesRangeNoI) {
        nCl2 *= (cl2Reg[i].end() - 1 - cl2Reg[i].begin());
    }

    // Compute the flux for the 0th order moments
    double f = _coefs(0, 0, 0, 0) * concentrations[_reactants[0]] *
        concentrations[_reactants[1]];
    for (auto i : speciesRangeNoI) {
        f += _coefs(i() + 1, 0, 0, 0) *
            concentrations[_reactantMomentIds[0][i()]] *
            concentrations[_reactants[1]];
    }
    for (auto j : speciesRangeNoI) {
        f += _coefs(0, j() + 1, 0, 0) * concentrations[_reactants[0]] *
            concentrations[_reactantMomentIds[1][j()]];
    }
    for (auto i : speciesRangeNoI) {
        for (auto j : speciesRangeNoI) {
            f += _coefs(i() + 1, j() + 1, 0, 0) *
                concentrations[_reactantMomentIds[0][i()]] *
                concentrations[_reactantMomentIds[1][j()]];
        }
    }
    f *= _rate(gridIndex);

    fluxes[_reactants[0]] -= f / (double) nCl1;
    fluxes[_reactants[1]] -= f / (double) nCl2;
    for (auto prodId : _products) {
        if (prodId == invalid) {
            continue;
        }

        auto prod = _network->getCluster(prodId);
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
        f = _coefs(0, 0, 1, k() + 1) * concentrations[_reactants[0]] *
            concentrations[_reactants[1]];
        for (auto i : speciesRangeNoI) {
            f += _coefs(i() + 1, 0, 1, k() + 1) *
                concentrations[_reactantMomentIds[0][i()]] *
                concentrations[_reactants[1]];
        }
        for (auto j : speciesRangeNoI) {
            f += _coefs(0, j() + 1, 1, k() + 1) *
                concentrations[_reactants[0]] *
                concentrations[_reactantMomentIds[1][j()]];
        }
        for (auto i : speciesRangeNoI) {
            for (auto j : speciesRangeNoI) {
                f += _coefs(i() + 1, j() + 1, 1, k() + 1) *
                    concentrations[_reactantMomentIds[0][i()]] *
                    concentrations[_reactantMomentIds[1][j()]];
            }
        }
        f *= _rate(gridIndex);
        fluxes[_reactantMomentIds[0][k()]] -= f / (double) nCl1;

        // For the second reactant
        f = _coefs(0, 0, 2, k() + 1) * concentrations[_reactants[0]] *
            concentrations[_reactants[1]];
        for (auto i : speciesRangeNoI) {
            f += _coefs(i() + 1, 0, 2, k() + 1) *
                concentrations[_reactantMomentIds[0][i()]] *
                concentrations[_reactants[1]];
        }
        for (auto j : speciesRangeNoI) {
            f += _coefs(0, j() + 1, 2, k() + 1) *
                concentrations[_reactants[0]] *
                concentrations[_reactantMomentIds[1][j()]];
        }
        for (auto i : speciesRangeNoI) {
            for (auto j : speciesRangeNoI) {
                f += _coefs(i() + 1, j() + 1, 2, k() + 1) *
                    concentrations[_reactantMomentIds[0][i()]] *
                    concentrations[_reactantMomentIds[1][j()]];
            }
        }
        f *= _rate(gridIndex);
        fluxes[_reactantMomentIds[1][k()]] -= f / (double) nCl2;

        // For the products
        for (auto prodId : _products) {
            if (prodId == invalid) {
                continue;
            }

            auto prod = _network->getCluster(prodId);
            const auto& prodReg = prod.getRegion();
            AmountType nProd = 1;
            for (auto i : speciesRangeNoI) {
                nProd *= (prodReg[i].end() - 1 - prodReg[i].begin());
            }

            // TODO generalize for 0, 1, or 2 products because here we are
            // taking the same coefs for all of them
            f = _coefs(0, 0, 0, k() + 1) * concentrations[_reactants[0]] *
                concentrations[_reactants[1]];
            for (auto i : speciesRangeNoI) {
                f += _coefs(i() + 1, 0, 0, k() + 1) *
                    concentrations[_reactantMomentIds[0][i()]] *
                    concentrations[_reactants[1]];
            }
            for (auto j : speciesRangeNoI) {
                f += _coefs(0, j() + 1, 0, k() + 1) *
                    concentrations[_reactants[0]] *
                    concentrations[_reactantMomentIds[1][j()]];
            }
            for (auto i : speciesRangeNoI) {
                for (auto j : speciesRangeNoI) {
                    f += _coefs(i() + 1, j() + 1, 0, k() + 1) *
                        concentrations[_reactantMomentIds[0][i()]] *
                        concentrations[_reactantMomentIds[1][j()]];
                }
            }
            f *= _rate(gridIndex);
            // TODO use the correct first index
            fluxes[_productMomentIds[0][k()]] -= f / (double) nProd;
        }
    }
}

template <typename TImpl>
template <typename TDerived>
inline void
ReactionNetwork<TImpl>::Reaction<TDerived>::dissociationFlux(
    ConcentrationsView concentrations, FluxesView fluxes, std::size_t gridIndex)
{
    using AmountType = typename NetworkType::AmountType;
    constexpr auto numSpeciesNoI = NetworkType::getNumberOfSpeciesNoI();
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

    // Compute the total number of elements in each cluster
    auto cl = _network->getCluster(_reactants[0]);
    const auto& clReg = cl.getRegion();
    AmountType nCl = 1;
    for (auto i : speciesRangeNoI) {
        nCl *= (clReg[i].end() - 1 - clReg[i].begin());
    }
    auto prod1 = _network->getCluster(_products[0]);
    const auto& prod1Reg = prod1.getRegion();
    AmountType nProd1 = 1;
    for (auto i : speciesRangeNoI) {
        nProd1 *= (prod1Reg[i].end() - 1 - prod1Reg[i].begin());
    }
    auto prod2 = _network->getCluster(_products[1]);
    const auto& prod2Reg = prod2.getRegion();
    AmountType nProd2 = 1;
    for (auto i : speciesRangeNoI) {
        nProd2 *= (prod2Reg[i].end() - 1 - prod2Reg[i].begin());
    }

    // Compute the flux for the 0th order moments
    double f = _coefs(0, 0, 0, 0) * concentrations[_reactants[0]];
    for (auto i : speciesRangeNoI) {
        f += _coefs(i() + 1, 0, 0, 0) *
            concentrations[_reactantMomentIds[0][i()]];
    }
    f *= _rate(gridIndex);
    fluxes[_reactants[0]] -= f / (double) nCl;
    fluxes[_products[0]] += f / (double) nProd1;
    fluxes[_products[1]] += f / (double) nProd2;

    // Take care of the first moments
    for (auto k : speciesRangeNoI) {
        // First for the reactant
        f = _coefs(0, 0, 0, k() + 1) * concentrations[_reactants[0]];
        for (auto i : speciesRangeNoI) {
            f += _coefs(i() + 1, 0, 0, k() + 1) *
                concentrations[_reactantMomentIds[0][i()]];
        }
        // TODO compute the prefactor related to the dispersion, it can be
        // moved to the coefs maybe
        double prefactor = 1.0;
        f *= _rate(gridIndex) * prefactor;
        fluxes[_reactantMomentIds[0][k()]] -= f / (double) nCl;

        // Now the first product
        f = _coefs(0, 0, 1, k() + 1) * concentrations[_reactants[0]];
        for (auto i : speciesRangeNoI) {
            f += _coefs(i() + 1, 0, 1, k() + 1) *
                concentrations[_reactantMomentIds[0][i()]];
        }
        // TODO compute the prefactor related to the dispersion, it can be
        // moved to the coefs maybe
        prefactor = 1.0;
        f *= _rate(gridIndex) * prefactor;
        fluxes[_productMomentIds[0][k()]] += f / (double) nProd1;

        // Finally the second product
        f = _coefs(0, 0, 2, k() + 1) * concentrations[_reactants[0]];
        for (auto i : speciesRangeNoI) {
            f += _coefs(i() + 1, 0, 2, k() + 1) *
                concentrations[_reactantMomentIds[0][i()]];
        }
        // TODO compute the prefactor related to the dispersion, it can be
        // moved to the coefs maybe
        prefactor = 1.0;
        f *= _rate(gridIndex) * prefactor;
        fluxes[_productMomentIds[1][k()]] += f / (double) nProd2;
    }
}
}
}
