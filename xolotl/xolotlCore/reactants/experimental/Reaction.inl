#pragma once

#include <algorithm>
#include <array>

#include <plsm/EnumIndexed.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TDerived>
template <typename TReactionNetwork>
Reaction<TDerived>::Reaction(TReactionNetwork& network, Type reactionType,
    std::size_t cluster0, std::size_t cluster1, std::size_t cluster2,
    std::size_t cluster3)
    :
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
    _rate(static_cast<TDerived*>(this)->computeRate(network))
{
    if (_type == Type::production) {
        //TODO: Vary third extent based number of valid products
        _coefs = Kokkos::View<double****>("Flux Coefficients", 5, 5, 3, 5);
        computeProductionCoefficients(network);
    }
    else {
        _coefs = Kokkos::View<double****>("Flux Coefficients", 5, 1, 3, 5);
        computeDissociationCoefficients(network);
    }
}

template <typename TDerived>
template <typename TCluster>
inline typename TCluster::NetworkType::AmountType
Reaction<TDerived>::computeOverlap(
    TCluster singleCl, TCluster pairCl1, TCluster pairCl2)
{
    using NetworkType = typename TCluster::NetworkType;
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

template <typename TDerived>
template <typename TReactionNetwork>
inline void
Reaction<TDerived>::computeProductionCoefficients(TReactionNetwork& network)
{
    using Species = typename TReactionNetwork::Species;
    using SpeciesSeq = typename TReactionNetwork::SpeciesSequence;
    using AmountType = typename TReactionNetwork::AmountType;

    // Find the overlap for this reaction
    constexpr auto numSpeciesNoI = TReactionNetwork::getNumberOfSpeciesNoI();
    constexpr auto speciesRangeNoI = TReactionNetwork::getSpeciesRangeNoI();

    auto cl1 = network.getCluster(_reactants[0]);
    auto cl2 = network.getCluster(_reactants[1]);
    auto prod = network.getCluster(_products[0]);

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

template <typename TDerived>
template <typename TReactionNetwork>
inline void
Reaction<TDerived>::computeDissociationCoefficients(TReactionNetwork& network)
{
    using Species = typename TReactionNetwork::Species;
    using SpeciesSeq = typename TReactionNetwork::SpeciesSequence;
    using AmountType = typename TReactionNetwork::AmountType;

    constexpr auto numSpeciesNoI = TReactionNetwork::getNumberOfSpeciesNoI();
    constexpr auto speciesRangeNoI = TReactionNetwork::getSpeciesRangeNoI();

    auto cl = network.getCluster(_reactants[0]);
    auto prod1 = network.getCluster(_products[0]);
    auto prod2 = network.getCluster(_products[1]);

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
}
}
