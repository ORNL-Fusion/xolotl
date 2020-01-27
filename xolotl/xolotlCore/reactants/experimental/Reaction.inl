#pragma once

#include <algorithm>
#include <array>

#include <plsm/EnumIndexed.h>

#include <Constants.h>
#include <MathUtils.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
Reaction<TNetwork, TDerived>::Reaction(detail::ReactionDataRef reactionData,
    ClusterDataRef clusterData, std::size_t reactionId, Type reactionType,
    std::size_t cluster0, std::size_t cluster1, std::size_t cluster2,
    std::size_t cluster3)
    :
    _clusterData(clusterData),
    _type(reactionType),
    _connectFn(
        _type == Type::production ? &Reaction::productionConnectivity :
            &Reaction::dissociationConnectivity),
    _fluxFn(
        _type == Type::production ? &Reaction::productionFlux :
            &Reaction::dissociationFlux),
    _partialsFn(
        _type == Type::production ? &Reaction::productionPartialDerivatives :
            &Reaction::dissociationPartialDerivatives),
    _reactants(
        _type == Type::production ? Kokkos::Array<std::size_t, 2>( {cluster0,
            cluster1}) : Kokkos::Array<std::size_t, 2>( {cluster0, invalid})),
    _products(_type == Type::production ? Kokkos::Array<std::size_t, 2>( {
        cluster2, cluster3}) : Kokkos::Array<std::size_t, 2>( {cluster1,
        cluster2})),
    _rate(reactionData.getRates(reactionId)),
    _coefs(reactionData.getCoefficients(reactionId)),
    _inverseMap(reactionData.inverseMap)
{
    for (std::size_t i : {0, 1}) {
        copyMomentIds(_reactants[i], _reactantMomentIds[i]);
        copyMomentIds(_products[i], _productMomentIds[i]);
    }

    if (_type == Type::production) {
        computeProductionCoefficients();
    }
    else {
        computeDissociationCoefficients();
    }

    updateRates();
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::updateRates()
{
    if (_type == Type::production) {
        computeProductionRates();
    }
    else {
        computeDissociationRates();
    }
}

template <typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<hasInterstitial<typename TRegion::EnumIndex>(),
    typename TRegion::ScalarType>
getISizeForOverlap(const TRegion& singleClReg, const TRegion& pairCl1Reg,
    const TRegion& pairCl2Reg)
{
    using Species = typename TRegion::EnumIndex;
    return pairCl1Reg[Species::I].begin() + pairCl2Reg[Species::I].begin() -
        singleClReg[Species::I].begin();
}

template <typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<!hasInterstitial<typename TRegion::EnumIndex>(),
    typename TRegion::ScalarType>
getISizeForOverlap(const TRegion& singleClReg, const TRegion& pairCl1Reg,
    const TRegion& pairCl2Reg)
{
    return 0;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
typename Reaction<TNetwork, TDerived>::AmountType
Reaction<TNetwork, TDerived>::computeOverlap(const Region& singleClReg,
    const Region& pairCl1Reg, const Region& pairCl2Reg)
{
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

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

        // Special case for I
        if (isVacancy(i)) {
            auto iSize = getISizeForOverlap(singleClReg, pairCl1Reg, pairCl2Reg);
            for (auto j : makeIntervalRange(pairCl1Reg[i])) {
                width +=
                    min(singleClReg[i].end() - 1, pairCl2Reg[i].end() - 1 + j - iSize)
                    - max(singleClReg[i].begin(), pairCl2Reg[i].begin() + j - iSize)
                    + 1;
            }
        }
        else {
            //TODO: Would be nice to loop on the cluster with the smaller tile
            for (auto j : makeIntervalRange(pairCl1Reg[i])) {
                width +=
                    min(singleClReg[i].end() - 1, pairCl2Reg[i].end() - 1 + j)
                    - max(singleClReg[i].begin(), pairCl2Reg[i].begin() + j)
                    + 1;
            }
        }
        nOverlap *= width;
    }

//    if (nOverlap <= 0) {
//        std::cout << pairCl1Reg[TNetwork::Traits::Species::He].begin() << ", " << pairCl1Reg[TNetwork::Traits::Species::D].begin()
//              << ", " << pairCl1Reg[TNetwork::Traits::Species::T].begin()
//            << ", " << pairCl1Reg[TNetwork::Traits::Species::V].begin() << ", " << pairCl1Reg[TNetwork::Traits::Species::I].begin() << std::endl;
//        std::cout << pairCl2Reg[TNetwork::Traits::Species::He].begin() << ", " << pairCl2Reg[TNetwork::Traits::Species::D].begin()
//              << ", " << pairCl2Reg[TNetwork::Traits::Species::T].begin()
//            << ", " << pairCl2Reg[TNetwork::Traits::Species::V].begin() << ", " << pairCl2Reg[TNetwork::Traits::Species::I].begin() << std::endl;
//        std::cout << "Prod: " << singleClReg[TNetwork::Traits::Species::He].begin() << ", " << singleClReg[TNetwork::Traits::Species::D].begin()
//              << ", " << singleClReg[TNetwork::Traits::Species::T].begin()
//            << ", " << singleClReg[TNetwork::Traits::Species::V].begin() << ", " << singleClReg[TNetwork::Traits::Species::I].begin() << std::endl;
//        std::cout << "Overlap: " << nOverlap << std::endl;
//    }
    assert(nOverlap > 0);

    return nOverlap;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::computeProductionCoefficients()
{
    // static
    const auto dummyRegion = Region(Composition{});

    // Find the overlap for this reaction
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

    const auto& cl1Reg = _clusterData.getCluster(_reactants[0]).getRegion();
    const auto& cl2Reg = _clusterData.getCluster(_reactants[1]).getRegion();
    const auto& prod1Reg = (_products[0] == invalid) ? dummyRegion :
        _clusterData.getCluster(_products[0]).getRegion();
    const auto& prod2Reg = (_products[1] == invalid) ? dummyRegion :
        _clusterData.getCluster(_products[1]).getRegion();
    const auto& cl1Disp = cl1Reg.dispersion();
    const auto& cl2Disp = cl2Reg.dispersion();
    
    // If there is no product the overlap is 1
    double nOverlap = 1.0;
    // General case
    if (_products[0] != invalid && _products[1] == invalid)
        nOverlap = static_cast<double>(computeOverlap(prod1Reg, cl1Reg, cl2Reg));
    // Special case with two products
    else if (_products[0] != invalid && _products[1] != invalid) {
        // Combine the regions
        auto ilist = Kokkos::Array<plsm::Interval<AmountType>, NetworkType::getNumberOfSpecies()>();
        for (auto i : NetworkType::getSpeciesRange()) {
            auto inter = plsm::Interval<AmountType> (prod1Reg[i].begin() + prod2Reg[i].begin(),
                    prod1Reg[i].end() + prod2Reg[i].end() - 1);
            ilist[i()] = inter;
        }
        auto prodReg = Region(ilist);
        nOverlap = static_cast<double>(computeOverlap(prodReg, cl1Reg, cl2Reg));
    }

    _coefs(0, 0, 0, 0) = nOverlap;
    for (auto i : speciesRangeNoI) {
        // First order sum on the first reactant
        for (double m : makeIntervalRange(prod2Reg[i]))
        for (double l : makeIntervalRange(cl2Reg[i])) {
            _coefs(i() + 1, 0, 0, 0) += firstOrderSum(
                max(prod1Reg[i].begin() + m - l, static_cast<double>(cl1Reg[i].begin())),
                min(prod1Reg[i].end() - 1 + m - l, static_cast<double>(cl1Reg[i].end() - 1)),
                static_cast<double>(cl1Reg[i].end() - 1 + cl1Reg[i].begin())
                    / 2.0);
        }

        // First order sum on the second reactant
        for (double m : makeIntervalRange(prod2Reg[i]))
        for (double l : makeIntervalRange(cl1Reg[i])) {
            _coefs(0, i() + 1, 0, 0) += firstOrderSum(
                max(prod1Reg[i].begin() + m - l, static_cast<double>(cl2Reg[i].begin())),
                min(prod1Reg[i].end() - 1 + m - l, static_cast<double>(cl2Reg[i].end() - 1)),
                static_cast<double>(cl2Reg[i].end() - 1 + cl2Reg[i].begin())
                    / 2.0);
        }

        // Loop on the potential products
        for (std::size_t p : {0,1}) {
            auto prodId = _products[p];
            if (prodId == invalid) {
                continue;
            }

            // Get the regions in the right order
            const auto& thisReg = (prodId == _products[0]) ? prod1Reg : prod2Reg;
            const auto& otherReg = (prodId == _products[0]) ? prod2Reg : prod1Reg;
            // Get the dispersion
            const auto& thisDispersion = thisReg.dispersion();

            // First order sum on the other product
            for (double m : makeIntervalRange(otherReg[i]))
            for (double l : makeIntervalRange(cl1Reg[i])) {
                _coefs(0, 0, p+2, i() + 1) += firstOrderSum( // p+2 because 0 and 1 are used for reactants
                    max(static_cast<double>(thisReg[i].begin()), cl2Reg[i].begin() + l - m),
                    min(static_cast<double>(thisReg[i].end() - 1), cl2Reg[i].end() - 1 + l - m),
                    static_cast<double>(thisReg[i].end() - 1 + thisReg[i].begin())
                        / 2.0);
            }
            _coefs(0, 0, p+2, i() + 1) /= thisDispersion[i()];

            // Products first moments
            for (auto k : speciesRangeNoI) {
                // Second order sum
                if (k == i) {
                    for (double m : makeIntervalRange(otherReg[i]))
                    for (double l : makeIntervalRange(cl2Reg[i])) {
                    _coefs(i() + 1, 0, p+2, k() + 1) += secondOrderOffsetSum(
                            max(thisReg[i].begin() + m - l, static_cast<double>(cl1Reg[i].begin())),
                            min(thisReg[i].end() - 1 + m - l, static_cast<double>(cl1Reg[i].end() - 1)),
                            static_cast<double>(cl1Reg[i].end() - 1 + cl1Reg[i].begin())
                                / 2.0,
                                static_cast<double>(thisReg[i].end() - 1 + thisReg[i].begin())
                                / 2.0, l - m);
                    }
                    _coefs(i() + 1, 0, p+2, k() + 1) /= thisDispersion[k()];

                    for (double m : makeIntervalRange(otherReg[i]))
                    for (double l : makeIntervalRange(cl1Reg[i])) {
                        _coefs(0, i() + 1, p+2, k() + 1) += secondOrderOffsetSum(
                            max(thisReg[i].begin() + m - l, static_cast<double>(cl2Reg[i].begin())),
                            min(thisReg[i].end() - 1 + m - l, static_cast<double>(cl2Reg[i].end() - 1)),
							static_cast<double>(cl2Reg[i].end() - 1 + cl2Reg[i].begin())
                                / 2.0,
								static_cast<double>(thisReg[i].end() - 1 + thisReg[i].begin())
                                / 2.0, l - m);
                    }
                    _coefs(0, i() + 1, p+2, k() + 1) /= thisDispersion[k()];

                }
                else {
                    _coefs(i() + 1, 0, p+2, k() + 1) = _coefs(i() + 1, 0, 0, 0)
                        * _coefs(0, 0, p+2, k() + 1) / nOverlap;

                    _coefs(0, i() + 1, p+2, k() + 1) = _coefs(0, i() + 1, 0, 0)
                        * _coefs(0, 0, p+2, k() + 1) / nOverlap;
                }
            }
        }
    }

    for (auto i : speciesRangeNoI) {
        // First reactant first moments
        for (auto k : speciesRangeNoI) {
            _coefs(0, 0, 0, k() + 1) = _coefs(k() + 1, 0, 0, 0) / cl1Disp[k()];

            if (k == i) {
                for (double m : makeIntervalRange(prod2Reg[i]))
                for (double l : makeIntervalRange(cl2Reg[i])) {
                    _coefs(i() + 1, 0, 0, k() + 1) += secondOrderSum(
                        max(prod1Reg[i].begin() + m - l, static_cast<double>(cl1Reg[i].begin())),
                        min(prod1Reg[i].end() - 1 + m - l, static_cast<double>(cl1Reg[i].end() - 1)),
						static_cast<double>(cl1Reg[i].end() - 1 + cl1Reg[i].begin())
                            / 2.0);
                }
                _coefs(i() + 1, 0, 0, k() + 1) /= cl1Disp[k()];
            }
            else {
                _coefs(i() + 1, 0, 0, k() + 1) = _coefs(i() + 1, 0, 0, 0)
                    * _coefs(k() + 1, 0, 0, 0) / (nOverlap * cl1Disp[k()]);
            }

            _coefs(0, i() + 1, 0, k() + 1) = _coefs(k() + 1, i() + 1, 0, 0) / cl1Disp[k()];
        }

        // Second reactant partial derivatives
        for (auto k : speciesRangeNoI) {
            _coefs(0, 0, 1, k() + 1) = _coefs(0, k() + 1, 0, 0) / cl2Disp[k()];

            if (k == i) {
                for (double m : makeIntervalRange(prod2Reg[i]))
                for (double l : makeIntervalRange(cl1Reg[i])) {
                    _coefs(0, i() + 1, 1, k() + 1) += secondOrderSum(
                        max(prod1Reg[i].begin() + m - l, static_cast<double>(cl2Reg[i].begin())),
                        min(prod1Reg[i].end() - 1 + m - l, static_cast<double>(cl2Reg[i].end() - 1)),
						static_cast<double>(cl2Reg[i].end() - 1 + cl2Reg[i].begin())
                            / 2.0);
                }
                _coefs(0, i() + 1, 1, k() + 1) /= cl2Disp[k()];
            }
            else {
                _coefs(0, i() + 1, 1, k() + 1) = _coefs(0, i() + 1, 0, 0)
                    * _coefs(0, k() + 1, 0, 0) / (nOverlap * cl2Disp[k()]);
            }

            _coefs(i() + 1, 0, 1, k() + 1) = _coefs(i() + 1, k() + 1, 0, 0) / cl2Disp[k()];
        }
    }

    // Now we loop over the 2 dimensions of the coefs to compute all
    // the possible sums over distances for the flux
    for (auto i : speciesRangeNoI) {
        for (auto j : speciesRangeNoI) {
            // Second order sum
            if (i == j) {
                for (double m : makeIntervalRange(prod2Reg[j]))
                for (double l : makeIntervalRange(cl1Reg[j])) {
                    _coefs(i() + 1, j() + 1, 0, 0) += (l
                        - static_cast<double>(cl1Reg[j].end() - 1 + cl1Reg[j].begin())
                            / 2.0)
                        * firstOrderSum(
                            max(prod1Reg[j].begin() + m - l, static_cast<double>(cl2Reg[j].begin())),
                            min(prod1Reg[j].end() - 1 + m - l,
                            		static_cast<double>(cl2Reg[j].end() - 1)),
									static_cast<double>(cl2Reg[j].end() - 1 + cl2Reg[j].begin())
                                / 2.0);
                }
            }
            else {
                _coefs(i() + 1, j() + 1, 0, 0) = _coefs(i() + 1, 0, 0, 0)
                    * _coefs(0, j() + 1, 0, 0) / nOverlap;
            }

            // Now we deal with the coefficients needed for the
            // first moments
            // Let's start with the products
            for (std::size_t p : {0,1}) {
                auto prodId = _products[p];
                if (prodId == invalid) {
                    continue;
                }

                // Get the regions in the right order
                const auto& thisReg = (prodId == _products[0]) ? prod1Reg : prod2Reg;
                const auto& otherReg = (prodId == _products[0]) ? prod2Reg : prod1Reg;
                // Get the dispersion
                const auto& thisDispersion = thisReg.dispersion();

                for (auto k : speciesRangeNoI) {
                    // Third order sum
                    if (i == j && j == k) {
                        for (double m : makeIntervalRange(otherReg[i]))
                        for (double l : makeIntervalRange(cl1Reg[i])) {
                            _coefs(i() + 1, j() + 1, p+2, k() + 1) += (l
                                - static_cast<double>(cl1Reg[i].end() - 1 + cl1Reg[i].begin())
                                    / 2.0)
                                * secondOrderOffsetSum(
                                    max(thisReg[i].begin() + m - l,
                                    		static_cast<double>(cl2Reg[i].begin())),
                                    min(thisReg[i].end() - 1 + m - l,
                                    		static_cast<double>(cl2Reg[i].end() - 1)),
											static_cast<double>(cl2Reg[i].end() - 1
                                        + cl2Reg[i].begin()) / 2.0,
										static_cast<double>(thisReg[i].end() - 1
                                        + thisReg[i].begin()) / 2.0, l - m);
                        }
                        _coefs(i() + 1, j() + 1, p+2, k() + 1) /= thisDispersion[k()];
                    }
                    else if (j == k) {
                        _coefs(i() + 1, j() + 1, p+2, k() + 1) = _coefs(i() + 1, 0,
                            0, 0) * _coefs(0, j() + 1, p+2, k() + 1) / nOverlap;
                    }
                    else if (i == k) {
                        _coefs(i() + 1, j() + 1, p+2, k() + 1) = _coefs(0, j() + 1,
                            0, 0) * _coefs(i() + 1, 0, p+2, k() + 1) / nOverlap;
                    }
                    else {
                        // TODO check this is the right formula, might be divided by nOverlap^2
                        _coefs(i() + 1, j() + 1, p+2, k() + 1) = _coefs(i() + 1, 0,
                            0, 0) * _coefs(0, j() + 1, 0, 0)
                            * _coefs(0, 0, p+2, k() + 1) / nOverlap;
                    }
                }
            }

            // Let's take care of the first reactant first moments
            for (auto k : speciesRangeNoI) {
                // Third order sum
                if (i == j && j == k) {
                    for (double m : makeIntervalRange(prod2Reg[i]))
                    for (double l : makeIntervalRange(cl1Reg[i])) {
                        _coefs(i() + 1, j() + 1, 0, k() + 1) += (l
                            - static_cast<double>(cl1Reg[i].end() - 1 + cl1Reg[i].begin())
                                / 2.0)
                            * (l
                                - static_cast<double>(cl1Reg[i].end() - 1
                                    + cl1Reg[i].begin()) / 2.0)
                            * firstOrderSum(
                                max(prod1Reg[i].begin() + m - l,
                                		static_cast<double>(cl2Reg[i].begin())),
                                min(prod1Reg[i].end() - 1 + m - l,
                                		static_cast<double>(cl2Reg[i].end() - 1)),
										static_cast<double>(cl2Reg[i].end() - 1
                                    + cl2Reg[i].begin()) / 2.0);
                    }
                    _coefs(i() + 1, j() + 1, 0, k() + 1) /= cl1Disp[k()];
                }
                else if (i == k) {
                    _coefs(i() + 1, j() + 1, 0, k() + 1) = _coefs(0, j() + 1,
                        0, 0) * _coefs(i() + 1, 0, 0, k() + 1) / nOverlap;
                }
                else if (j == k) {
                    _coefs(i() + 1, j() + 1, 0, k() + 1) = _coefs(i() + 1, 0,
                        0, 0) * _coefs(0, j() + 1, 0, k() + 1) / nOverlap;
                }
                else {
                    // TODO check this is the right formula, might be divided by nOverlap^2
                    _coefs(i() + 1, j() + 1, 0, k() + 1) = _coefs(i() + 1, 0,
                        0, 0) * _coefs(0, j() + 1, 0, 0)
                        * _coefs(k() + 1, 0, 0, 0) / (nOverlap * cl1Disp[k()]);
                }
            }

            // Let's take care of the second reactant partial derivatives
            for (auto k : speciesRangeNoI) {
                // Third order sum
                if (i == j && j == k) {
                    for (double m : makeIntervalRange(prod2Reg[i]))
                    for (double l : makeIntervalRange(cl2Reg[i])) {
                        _coefs(i() + 1, j() + 1, 1, k() + 1) += (l
                            - static_cast<double>(cl2Reg[i].end() - 1 + cl2Reg[i].begin())
                                / 2.0)
                            * (l
                                - static_cast<double>(cl2Reg[i].end() - 1
                                    + cl2Reg[i].begin()) / 2.0)
                            * firstOrderSum(
                                max(prod1Reg[i].begin() + m - l,
                                		static_cast<double>(cl1Reg[i].begin())),
                                min(prod1Reg[i].end() - 1 + m - l,
                                		static_cast<double>(cl1Reg[i].end() - 1)),
                                (double) (cl1Reg[i].end() - 1
                                    + cl1Reg[i].begin()) / 2.0);
                    }
                    _coefs(i() + 1, j() + 1, 1, k() + 1) /= cl2Disp[k()];
                }
                else if (i == k) {
                    _coefs(i() + 1, j() + 1, 1, k() + 1) = _coefs(0, j() + 1,
                        0, 0) * _coefs(i() + 1, 0, 1, k() + 1) / nOverlap;
                }
                else if (j == k) {
                    _coefs(i() + 1, j() + 1, 1, k() + 1) = _coefs(i() + 1, 0,
                        0, 0) * _coefs(0, j() + 1, 1, k() + 1) / nOverlap;
                }
                else {
                    // TODO check this is the right formula, might be divided by nOverlap^2
                    _coefs(i() + 1, j() + 1, 1, k() + 1) = _coefs(i() + 1, 0,
                        0, 0) * _coefs(0, j() + 1, 0, 0)
                        * _coefs(0, k() + 1, 0, 0) / (nOverlap * cl2Disp[k()]);
                }
            }
        }
    }
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::computeDissociationCoefficients()
{
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

    auto clReg = _clusterData.getCluster(_reactants[0]).getRegion();
    auto prod1Reg = _clusterData.getCluster(_products[0]).getRegion();
    auto prod2Reg = _clusterData.getCluster(_products[1]).getRegion();
    const auto& clDisp = clReg.dispersion();
    const auto& prod1Disp = prod1Reg.dispersion();
    const auto& prod2Disp = prod2Reg.dispersion();

    auto nOverlap =
        static_cast<double>(computeOverlap(clReg, prod1Reg, prod2Reg));

    // The first coefficient is simply the overlap because it is the sum over 1
    _coefs(0, 0, 0, 0) = nOverlap;
    for (auto i : speciesRangeNoI) {
        // First order sum
        for (double l : makeIntervalRange(prod1Reg[i])) {
            _coefs(i() + 1, 0, 0, 0) += firstOrderSum(
                max(static_cast<double>(clReg[i].begin()), prod2Reg[i].begin() + l),
                min(static_cast<double>(clReg[i].end() - 1), prod2Reg[i].end() - 1 + l),
				static_cast<double>(clReg[i].end() - 1 + clReg[i].begin()) / 2.0);
        }
    }

    // First moments
    for (auto k : speciesRangeNoI) {
        // Reactant
        _coefs(0, 0, 0, k() + 1) = _coefs(k() + 1, 0, 0, 0) / clDisp[k()];

        // First product
        for (double l : makeIntervalRange(prod2Reg[k])) {
            _coefs(0, 0, 1, k() + 1) += firstOrderSum(
                max(clReg[k].begin() - l, static_cast<double>(prod1Reg[k].begin())),
                min(clReg[k].end() - 1 - l, static_cast<double>(prod1Reg[k].end() - 1)),
				static_cast<double>(prod1Reg[k].end() - 1 + prod1Reg[k].begin()) / 2.0);
        }
        _coefs(0, 0, 1, k() + 1) /= prod1Disp[k()];

        // Second product
        for (double l : makeIntervalRange(prod1Reg[k])) {
            _coefs(0, 0, 2, k() + 1) += firstOrderSum(
                max(clReg[k].begin() - l, static_cast<double>(prod2Reg[k].begin())),
                min(clReg[k].end() - 1 - l, static_cast<double>(prod2Reg[k].end() - 1)),
				static_cast<double>(prod2Reg[k].end() - 1 + prod2Reg[k].begin()) / 2.0);
        }
        _coefs(0, 0, 2, k() + 1) /= prod2Disp[k()];
    }

    // Now we loop over the 1 dimension of the coefs to compute all the
    // possible sums over distances for the flux
    for (auto i : speciesRangeNoI) {
        // Now we deal with the coefficients needed for the partial derivatives
        // Starting with the reactant
        for (auto k : speciesRangeNoI) {
            // Second order sum
            if (k == i) {
                for (double l : makeIntervalRange(prod1Reg[i])) {
                    _coefs(i() + 1, 0, 0, k() + 1) += secondOrderSum(
                        max(static_cast<double>(clReg[i].begin()), prod2Reg[i].begin() + l),
                        min(static_cast<double>(clReg[i].end() - 1), prod2Reg[i].end() - 1 + l),
						static_cast<double>(clReg[i].end() - 1 + clReg[i].begin()) / 2.0);
                }
                _coefs(i() + 1, 0, 0, k() + 1) /= clDisp[k()];
            }
            else {
                _coefs(i() + 1, 0, 0, k() + 1) = _coefs(i() + 1, 0, 0, 0)
                    * _coefs(k() + 1, 0, 0, 0) / (nOverlap * clDisp[k()]);
            }
        }

        // First moments for the first product
        for (auto k : speciesRangeNoI) {
            // Second order sum
            if (k == i) {
                for (double l : makeIntervalRange(prod2Reg[i])) {
                    _coefs(i() + 1, 0, 1, k() + 1) += secondOrderOffsetSum(
                        max(static_cast<double>(clReg[i].begin()), prod1Reg[i].begin() + l),
                        min(static_cast<double>(clReg[i].end() - 1), prod1Reg[i].end() - 1 + l),
						static_cast<double>(clReg[i].end() - 1 + clReg[i].begin()) / 2.0,
						static_cast<double>(prod1Reg[i].end() - 1 + prod1Reg[i].begin())
                            / 2.0, -l);
                }
                _coefs(i() + 1, 0, 1, k() + 1) /= prod1Disp[k()];
            }
            else {
                _coefs(i() + 1, 0, 1, k() + 1) = _coefs(i() + 1, 0, 0, 0)
                    * _coefs(0, 0, 1, k() + 1) / nOverlap;
            }
        }

        // First moments for the second product
        for (auto k : speciesRangeNoI) {
            // Second order sum
            if (k == i) {
                for (double l : makeIntervalRange(prod1Reg[i])) {
                    _coefs(i() + 1, 0, 2, k() + 1) += secondOrderOffsetSum(
                        max(static_cast<double>(clReg[i].begin()), prod2Reg[i].begin() + l),
                        min(static_cast<double>(clReg[i].end() - 1), prod2Reg[i].end() - 1 + l),
						static_cast<double>(clReg[i].end() - 1 + clReg[i].begin()) / 2.0,
						static_cast<double>(prod2Reg[i].end() - 1 + prod2Reg[i].begin())
                            / 2.0, -l);
                }
                _coefs(i() + 1, 0, 2, k() + 1) /= prod2Disp[k()];
            }
            else {
                _coefs(i() + 1, 0, 2, k() + 1) = _coefs(i() + 1, 0, 0, 0)
                    * _coefs(0, 0, 2, k() + 1) / nOverlap;
            }
        }
    }
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
Reaction<TNetwork, TDerived>::computeProductionRate(std::size_t gridIndex)
{
    auto cl0 = _clusterData.getCluster(_reactants[0]);
    auto cl1 = _clusterData.getCluster(_reactants[1]);

    double r0 = cl0.getReactionRadius();
    double r1 = cl1.getReactionRadius();

    double dc0 = cl0.getDiffusionCoefficient(gridIndex);
    double dc1 = cl1.getDiffusionCoefficient(gridIndex);

    constexpr double pi = ::xolotlCore::pi;

    double kPlus = 4.0 * pi * (r0 + r1) * (dc0 + dc1);

    return kPlus;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
Reaction<TNetwork, TDerived>::computeDissociationRate(std::size_t gridIndex)
{
    double omega = _clusterData.getAtomicVolume();
    double T = _clusterData.temperature(gridIndex);

    // TODO: computeProductionRate should use products and not reactants
    auto cl0 = _clusterData.getCluster(_products[0]);
    auto cl1 = _clusterData.getCluster(_products[1]);

    double r0 = cl0.getReactionRadius();
    double r1 = cl1.getReactionRadius();

    double dc0 = cl0.getDiffusionCoefficient(gridIndex);
    double dc1 = cl1.getDiffusionCoefficient(gridIndex);

    constexpr double pi = ::xolotlCore::pi;

    double kPlus = 4.0 * pi * (r0 + r1) * (dc0 + dc1);
    double E_b = asDerived()->computeBindingEnergy();

    constexpr double k_B = ::xolotlCore::kBoltzmann;

    double kMinus = (1.0 / omega) * kPlus * std::exp(-E_b / (k_B * T));

    return kMinus;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::productionConnectivity(
    ConnectivityView connectivity)
{
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
    // Get the total number of elements in each cluster
    auto cl1 = _clusterData.getCluster(_reactants[0]);
    const auto& cl1Reg = cl1.getRegion();
    const bool cl1IsSimplex = cl1Reg.isSimplex();
    auto cl2 = _clusterData.getCluster(_reactants[1]);
    const auto& cl2Reg = cl2.getRegion();
    const bool cl2IsSimplex = cl2Reg.isSimplex();
    // Each reactant connects with all the reactants
    // Reactant 1 with reactant 1
    addConnectivity(_reactants[0], _reactants[0], connectivity);
    if (!cl1IsSimplex) {
        for (auto i : speciesRangeNoI) {
            addConnectivity(_reactants[0], _reactantMomentIds[0][i()], connectivity);
            addConnectivity(_reactantMomentIds[0][i()], _reactants[0], connectivity);
            for (auto j : speciesRangeNoI) {
                addConnectivity(_reactantMomentIds[0][i()], _reactantMomentIds[0][j()], connectivity);
            }
        }
    }
    // Reactant 2 with reactant 1
    addConnectivity(_reactants[1], _reactants[0], connectivity);
    if (!cl1IsSimplex) {
        for (auto i : speciesRangeNoI) {
            addConnectivity(_reactants[1], _reactantMomentIds[0][i()], connectivity);
        }
    }
    if (!cl2IsSimplex) {
        for (auto i : speciesRangeNoI) {
            addConnectivity(_reactantMomentIds[1][i()], _reactants[0], connectivity);
        }
    }
    if (!cl1IsSimplex && !cl2IsSimplex) {
        for (auto i : speciesRangeNoI) {
            for (auto j : speciesRangeNoI) {
                addConnectivity(_reactantMomentIds[1][i()], _reactantMomentIds[0][j()], connectivity);
            }
        }
    }
    // Reactant 1 with reactant 2
    addConnectivity(_reactants[0], _reactants[1], connectivity);
    if (!cl2IsSimplex) {
        for (auto i : speciesRangeNoI) {
            addConnectivity(_reactants[0], _reactantMomentIds[1][i()], connectivity);
        }
    }
    if (!cl1IsSimplex) {
        for (auto i : speciesRangeNoI) {
            addConnectivity(_reactantMomentIds[0][i()], _reactants[1], connectivity);
        }
    }
    if (!cl1IsSimplex && !cl2IsSimplex) {
        for (auto i : speciesRangeNoI) {
            for (auto j : speciesRangeNoI) {
                addConnectivity(_reactantMomentIds[0][i()], _reactantMomentIds[1][j()], connectivity);
            }
        }
    }
    // Reactant 2 with reactant 2
    addConnectivity(_reactants[1], _reactants[1], connectivity);
    if (!cl2IsSimplex) {
        for (auto i : speciesRangeNoI) {
            addConnectivity(_reactants[1], _reactantMomentIds[1][i()], connectivity);
            addConnectivity(_reactantMomentIds[1][i()], _reactants[1], connectivity);
            for (auto j : speciesRangeNoI) {
                addConnectivity(_reactantMomentIds[1][i()], _reactantMomentIds[1][j()], connectivity);
            }
        }
    }
    // Each product connects with all the reactants
    for (std::size_t p : {0,1}) {
        auto prodId = _products[p];
        if (prodId == invalid) {
            continue;
        }
        auto prod = _clusterData.getCluster(prodId);
        const auto& prodReg = prod.getRegion();
        const bool prodIsSimplex = prodReg.isSimplex();

        // With reactant 1
        addConnectivity(prodId, _reactants[0], connectivity);
        if (!cl1IsSimplex) {
            for (auto i : speciesRangeNoI) {
                addConnectivity(prodId, _reactantMomentIds[0][i()], connectivity);
            }
        }
        if (!prodIsSimplex) {
            for (auto i : speciesRangeNoI) {
                addConnectivity(_productMomentIds[p][i()], _reactants[0], connectivity);
            }
        }
        if (!cl1IsSimplex && !prodIsSimplex) {
            for (auto i : speciesRangeNoI) {
                for (auto j : speciesRangeNoI) {
                    addConnectivity(_productMomentIds[p][i()], _reactantMomentIds[0][j()], connectivity);
                }
            }
        }
        // With reactant 2
        addConnectivity(prodId, _reactants[1], connectivity);
        if (!cl2IsSimplex) {
            for (auto i : speciesRangeNoI) {
                addConnectivity(prodId, _reactantMomentIds[1][i()], connectivity);
            }
        }
        if (!prodIsSimplex) {
            for (auto i : speciesRangeNoI) {
                addConnectivity(_productMomentIds[p][i()], _reactants[1], connectivity);
            }
        }
        if (!cl2IsSimplex && !prodIsSimplex) {
            for (auto i : speciesRangeNoI) {
                for (auto j : speciesRangeNoI) {
                    addConnectivity(_productMomentIds[p][i()], _reactantMomentIds[1][j()], connectivity);
                }
            }
        }
    }
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::dissociationConnectivity(
    ConnectivityView connectivity)
{
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

    // Get the total number of elements in each cluster
    auto cl = _clusterData.getCluster(_reactants[0]);
    const auto& clReg = cl.getRegion();
    const bool clIsSimplex = clReg.isSimplex();
    auto prod1 = _clusterData.getCluster(_products[0]);
    const auto& prod1Reg = prod1.getRegion();
    const bool prod1IsSimplex = prod1Reg.isSimplex();
    auto prod2 = _clusterData.getCluster(_products[1]);
    const auto& prod2Reg = prod2.getRegion();
    const bool prod2IsSimplex = prod2Reg.isSimplex();

    // The reactant connects with the reactant
    addConnectivity(_reactants[0], _reactants[0], connectivity);
    if (!clIsSimplex) {
        for (auto i : speciesRangeNoI) {
            addConnectivity(_reactants[0], _reactantMomentIds[0][i()], connectivity);
            addConnectivity(_reactantMomentIds[0][i()], _reactants[0], connectivity);
            for (auto j : speciesRangeNoI) {
                addConnectivity(_reactantMomentIds[0][i()], _reactantMomentIds[0][j()], connectivity);
            }
        }
    }
    // Each product connects with  the reactant
    // Product 1 with reactant
    addConnectivity(_products[0], _reactants[0], connectivity);
    if (!clIsSimplex) {
        for (auto i : speciesRangeNoI) {
            addConnectivity(_products[0], _reactantMomentIds[0][i()], connectivity);
        }
    }
    if (!prod1IsSimplex) {
        for (auto i : speciesRangeNoI) {
            addConnectivity(_productMomentIds[0][i()], _reactants[0], connectivity);
        }
    }
    if (!clIsSimplex && !prod1IsSimplex) {
        for (auto i : speciesRangeNoI) {
            for (auto j : speciesRangeNoI) {
                addConnectivity(_productMomentIds[0][i()], _reactantMomentIds[0][j()], connectivity);
            }
        }
    }
    // Product 2 with reactant
    addConnectivity(_products[1], _reactants[0], connectivity);
    if (!clIsSimplex) {
        for (auto i : speciesRangeNoI) {
            addConnectivity(_products[1], _reactantMomentIds[0][i()], connectivity);
        }
    }
    if (!prod2IsSimplex) {
        for (auto i : speciesRangeNoI) {
            addConnectivity(_productMomentIds[1][i()], _reactants[0], connectivity);
        }
    }
    if (!clIsSimplex && !prod2IsSimplex) {
        for (auto i : speciesRangeNoI) {
            for (auto j : speciesRangeNoI) {
                addConnectivity(_productMomentIds[1][i()], _reactantMomentIds[0][j()], connectivity);
            }
        }
    }
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::productionFlux(ConcentrationsView concentrations,
    FluxesView fluxes, std::size_t gridIndex)
{
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

    // Compute the total number of elements in each cluster
    auto cl1 = _clusterData.getCluster(_reactants[0]);
    const auto& cl1Reg = cl1.getRegion();
    AmountType volCl1 = cl1Reg.volume();
    auto cl2 = _clusterData.getCluster(_reactants[1]);
    const auto& cl2Reg = cl2.getRegion();
    AmountType volCl2 = cl2Reg.volume();
    
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

    fluxes[_reactants[0]] -= f / (double) volCl1;
    fluxes[_reactants[1]] -= f / (double) volCl2;
    for (auto prodId : _products) {
        if (prodId == invalid) {
            continue;
        }

        auto prod = _clusterData.getCluster(prodId);
        const auto& prodReg = prod.getRegion();
        AmountType volProd = prodReg.volume();
        fluxes[prodId] += f / (double) volProd;
    }

    // Take care of the first moments
    for (auto k : speciesRangeNoI) {
        // First for the first reactant
    	if (volCl1 > 1) {
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
        fluxes[_reactantMomentIds[0][k()]] -= f / (double) volCl1;
    	}

        // For the second reactant
    	if (volCl2 > 1) {
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
        fluxes[_reactantMomentIds[1][k()]] -= f / (double) volCl2;
    	}

        // For the products
        for (std::size_t p : {0,1}) {
            auto prodId = _products[p];
            if (prodId == invalid) {
                continue;
            }

            auto prod = _clusterData.getCluster(prodId);
            const auto& prodReg = prod.getRegion();
            AmountType volProd = prodReg.volume();

            if (volProd > 1) {
            f = _coefs(0, 0, p+2, k() + 1) * concentrations[_reactants[0]] *
                concentrations[_reactants[1]];
            for (auto i : speciesRangeNoI) {
                f += _coefs(i() + 1, 0, p+2, k() + 1) *
                    concentrations[_reactantMomentIds[0][i()]] *
                    concentrations[_reactants[1]];
            }
            for (auto j : speciesRangeNoI) {
                f += _coefs(0, j() + 1, p+2, k() + 1) *
                    concentrations[_reactants[0]] *
                    concentrations[_reactantMomentIds[1][j()]];
            }
            for (auto i : speciesRangeNoI) {
                for (auto j : speciesRangeNoI) {
                    f += _coefs(i() + 1, j() + 1, p+2, k() + 1) *
                        concentrations[_reactantMomentIds[0][i()]] *
                        concentrations[_reactantMomentIds[1][j()]];
                }
            }
            f *= _rate(gridIndex);
            fluxes[_productMomentIds[p][k()]] += f / (double) volProd;
            }
        }
    }
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::dissociationFlux(
    ConcentrationsView concentrations, FluxesView fluxes, std::size_t gridIndex)
{
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

    // Compute the total number of elements in each cluster
    auto cl = _clusterData.getCluster(_reactants[0]);
    const auto& clReg = cl.getRegion();
    AmountType volCl = clReg.volume();
    auto prod1 = _clusterData.getCluster(_products[0]);
    const auto& prod1Reg = prod1.getRegion();
    AmountType volProd1 = prod1Reg.volume();
    auto prod2 = _clusterData.getCluster(_products[1]);
    const auto& prod2Reg = prod2.getRegion();
    AmountType volProd2 = prod2Reg.volume();

    // Compute the flux for the 0th order moments
    double f = _coefs(0, 0, 0, 0) * concentrations[_reactants[0]];
    for (auto i : speciesRangeNoI) {
        f += _coefs(i() + 1, 0, 0, 0) *
            concentrations[_reactantMomentIds[0][i()]];
    }
    f *= _rate(gridIndex);
    fluxes[_reactants[0]] -= f / (double) volCl;
    fluxes[_products[0]] += f / (double) volProd1;
    fluxes[_products[1]] += f / (double) volProd2;

    // Take care of the first moments
    for (auto k : speciesRangeNoI) {
        // First for the reactant
    	if (volCl > 1) {
        f = _coefs(0, 0, 0, k() + 1) * concentrations[_reactants[0]];
        for (auto i : speciesRangeNoI) {
            f += _coefs(i() + 1, 0, 0, k() + 1) *
                concentrations[_reactantMomentIds[0][i()]];
        }
        f *= _rate(gridIndex);
        fluxes[_reactantMomentIds[0][k()]] -= f / (double) volCl;
    	}

        // Now the first product
    	if (volProd1 > 1) {
        f = _coefs(0, 0, 1, k() + 1) * concentrations[_reactants[0]];
        for (auto i : speciesRangeNoI) {
            f += _coefs(i() + 1, 0, 1, k() + 1) *
                concentrations[_reactantMomentIds[0][i()]];
        }
        f *= _rate(gridIndex);
        fluxes[_productMomentIds[0][k()]] += f / (double) volProd1;
    	}

        // Finally the second product
    	if (volProd2 > 1) {
        f = _coefs(0, 0, 2, k() + 1) * concentrations[_reactants[0]];
        for (auto i : speciesRangeNoI) {
            f += _coefs(i() + 1, 0, 2, k() + 1) *
                concentrations[_reactantMomentIds[0][i()]];
        }
        f *= _rate(gridIndex);
        fluxes[_productMomentIds[1][k()]] += f / (double) volProd2;
    	}
    }
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::productionPartialDerivatives(
    ConcentrationsView concentrations,
    Kokkos::View<double*> values, std::size_t gridIndex)
{
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
    int nProd = 0;
    for (auto prodId : _products) {
        if (prodId != invalid) {
            ++nProd;
        }
    }

    // Compute the total number of elements in each cluster
    auto cl1 = _clusterData.getCluster(_reactants[0]);
    const auto& cl1Reg = cl1.getRegion();
    AmountType volCl1 = cl1Reg.volume();
    auto cl2 = _clusterData.getCluster(_reactants[1]);
    const auto& cl2Reg = cl2.getRegion();
    AmountType volCl2 = cl2Reg.volume();

    // Compute the partials for the 0th order moments
    // Compute the values (d / dL_0^A)
    double temp = _coefs(0, 0, 0, 0) * concentrations[_reactants[1]];
    if (volCl2 > 1) {
        for (auto i : speciesRangeNoI) {
            temp += _coefs(0, i() + 1, 0, 0) *
                concentrations[_reactantMomentIds[1][i()]];
        }
    }
    // First for the first reactant
    values(_inverseMap(_reactants[0], _reactants[0])) -= _rate(gridIndex) * temp / (double) volCl1;
    // Second reactant
    values(_inverseMap(_reactants[1], _reactants[0])) -= _rate(gridIndex) * temp / (double) volCl2;
    // For the products
    for (std::size_t p : {0,1}) {
        auto prodId = _products[p];
        if (prodId == invalid) {
            continue;
        }
        auto prod = _clusterData.getCluster(prodId);
        const auto& prodReg = prod.getRegion();
        AmountType volProd = prodReg.volume();
        values(_inverseMap(prodId, _reactants[0])) += _rate(gridIndex) * temp / (double) volProd;
    }
    
    // Compute the values (d / dL_0^B)
    temp = _coefs(0, 0, 0, 0) * concentrations[_reactants[0]];
    if (volCl1 > 1) {
        for (auto i : speciesRangeNoI) {
            temp += _coefs(i() + 1, 0, 0, 0) *
                concentrations[_reactantMomentIds[0][i()]];
        }
    }
    // First for the first reactant
    values(_inverseMap(_reactants[0], _reactants[1])) -= _rate(gridIndex) * temp / (double) volCl1;
    // Second reactant
    values(_inverseMap(_reactants[1], _reactants[1])) -= _rate(gridIndex) * temp / (double) volCl2;
    // For the products
    for (std::size_t p : {0,1}) {
        auto prodId = _products[p];
        if (prodId == invalid) {
            continue;
        }
        auto prod = _clusterData.getCluster(prodId);
        const auto& prodReg = prod.getRegion();
        AmountType volProd = prodReg.volume();
        values(_inverseMap(prodId, _reactants[1])) += _rate(gridIndex) * temp / (double) volProd;
    }
    
    // (d / dL_1^A)
    if (volCl1 > 1) {
    for (auto i : speciesRangeNoI) {
        temp = _coefs(i() + 1, 0, 0, 0) * concentrations[_reactants[1]];
        if (volCl2 > 1) {
        for (auto j : speciesRangeNoI) {
            temp += _coefs(i() + 1, j() + 1, 0, 0) * concentrations[_reactantMomentIds[1][j()]];
        }
        }
        // First reactant
        values(_inverseMap(_reactants[0], _reactantMomentIds[0][i()])) -= _rate(gridIndex) * temp / (double) volCl1;
        // second reactant
        values(_inverseMap(_reactants[1], _reactantMomentIds[0][i()])) -= _rate(gridIndex) * temp / (double) volCl2;
        // For the products
        for (std::size_t p : {0,1}) {
            auto prodId = _products[p];
            if (prodId == invalid) {
                continue;
            }
            auto prod = _clusterData.getCluster(prodId);
            const auto& prodReg = prod.getRegion();
            AmountType volProd = prodReg.volume();
            values(_inverseMap(prodId, _reactantMomentIds[0][i()])) += _rate(gridIndex) * temp / (double) volProd;
        }
    }
    }

    // (d / dL_1^B)
    if (volCl2 > 1) {
    for (auto i : speciesRangeNoI) {
        temp = _coefs(0, i() + 1, 0, 0) * concentrations[_reactants[0]];
        if (volCl1 > 1) {
        for (auto j : speciesRangeNoI) {
            temp += _coefs(j() + 1, i() + 1, 0, 0) * concentrations[_reactantMomentIds[0][j()]];
        }
        }
        values(_inverseMap(_reactants[0], _reactantMomentIds[1][i()])) -= _rate(gridIndex) * temp / (double) volCl1;
        values(_inverseMap(_reactants[1], _reactantMomentIds[1][i()])) -= _rate(gridIndex) * temp / (double) volCl2;
        for (std::size_t p : {0,1}) {
            auto prodId = _products[p];
            if (prodId == invalid) {
                continue;
            }
            auto prod = _clusterData.getCluster(prodId);
            const auto& prodReg = prod.getRegion();
            AmountType volProd = prodReg.volume();
            values(_inverseMap(prodId, _reactantMomentIds[1][i()])) += _rate(gridIndex) * temp / (double) volProd;
        }
    }
    }

    // Take care of the first moments
    if (volCl1 > 1) {
    for (auto k : speciesRangeNoI) {
        // First for the first reactant
        // (d / dL_0^A)
        temp = _coefs(0, 0, 0, k() + 1) * concentrations[_reactants[1]];
        if (volCl2 > 1) {
        for (auto j : speciesRangeNoI) {
            temp += _coefs(0, j() + 1, 0, k() + 1) * concentrations[_reactantMomentIds[1][j()]];
        }
        }
        values(_inverseMap(_reactantMomentIds[0][k()], _reactants[0])) -= _rate(gridIndex) * temp / (double) volCl1;
        // (d / dL_1^A)
        for (auto i : speciesRangeNoI) {
            temp = _coefs(i() + 1, 0, 0, k() + 1) * concentrations[_reactants[1]];
            if (volCl2 > 1) {
            for (auto j : speciesRangeNoI) {
                temp += _coefs(i() + 1, j() + 1, 0, k() + 1) * concentrations[_reactantMomentIds[1][j()]];
            }
            }
            values(_inverseMap(_reactantMomentIds[0][k()], _reactantMomentIds[0][i()]))
            -= _rate(gridIndex) * temp / (double) volCl1;
        }
        // (d / dL_0^B)
        temp = _coefs(0, 0, 0, k() + 1) * concentrations[_reactants[0]];
        for (auto j : speciesRangeNoI) {
            temp += _coefs(j() + 1, 0, 0, k() + 1) * concentrations[_reactantMomentIds[0][j()]];
        }
        values(_inverseMap(_reactantMomentIds[0][k()], _reactants[1]))
        -= _rate(gridIndex) * temp / (double) volCl1;
        // (d / dL_1^B)
        for (auto i : speciesRangeNoI) {
            temp = _coefs(0, i() + 1, 0, k() + 1) * concentrations[_reactants[0]];
            for (auto j : speciesRangeNoI) {
                temp += _coefs(j() + 1, i() + 1, 0, k() + 1) * concentrations[_reactantMomentIds[0][j()]];
            }
            values(_inverseMap(_reactantMomentIds[0][k()], _reactantMomentIds[1][i()]))
            -= _rate(gridIndex) * temp / (double) volCl1;
        }
    }
    }

    // Take care of the first moments
    if (volCl2 > 1) {
    for (auto k : speciesRangeNoI) {
        // First for the second reactant
        // (d / dL_0^A)
        temp = _coefs(0, 0, 1, k() + 1) * concentrations[_reactants[1]];
        for (auto j : speciesRangeNoI) {
            temp += _coefs(0, j() + 1, 1, k() + 1) * concentrations[_reactantMomentIds[1][j()]];
        }
        values(_inverseMap(_reactantMomentIds[1][k()], _reactants[0]))
        -= _rate(gridIndex) * temp / (double) volCl2;
        // (d / dL_1^A)
        for (auto i : speciesRangeNoI) {
            temp = _coefs(i() + 1, 0, 1, k() + 1) * concentrations[_reactants[1]];
            for (auto j : speciesRangeNoI) {
                temp += _coefs(i() + 1, j() + 1, 1, k() + 1) * concentrations[_reactantMomentIds[1][j()]];
            }
            values(_inverseMap(_reactantMomentIds[1][k()], _reactantMomentIds[0][i()]))
            -= _rate(gridIndex) * temp / (double) volCl2;
        }
        // (d / dL_0^B)
        temp = _coefs(0, 0, 1, k() + 1) * concentrations[_reactants[0]];
        if (volCl1 > 1) {
        for (auto j : speciesRangeNoI) {
            temp += _coefs(j() + 1, 0, 1, k() + 1) * concentrations[_reactantMomentIds[0][j()]];
        }
        }
        values(_inverseMap(_reactantMomentIds[1][k()], _reactants[1]))
        -= _rate(gridIndex) * temp / (double) volCl2;
        // (d / dL_1^B)
        for (auto i : speciesRangeNoI) {
            temp = _coefs(0, i() + 1, 1, k() + 1) * concentrations[_reactants[0]];
            if (volCl1 > 1) {
            for (auto j : speciesRangeNoI) {
                temp += _coefs(j() + 1, i() + 1, 1, k() + 1) * concentrations[_reactantMomentIds[0][j()]];
            }
            }
            values(_inverseMap(_reactantMomentIds[1][k()], _reactantMomentIds[1][i()]))
            -= _rate(gridIndex) * temp / (double) volCl2;
        }
    }
    }

    // Loop on the products
    for (std::size_t p : {0,1}) {
        auto prodId = _products[p];
        if (prodId == invalid) {
            continue;
        }

        auto prod = _clusterData.getCluster(prodId);
        const auto& prodReg = prod.getRegion();
        AmountType volProd = prodReg.volume();

        // Take care of the first moments
        if (volProd > 1) {
        for (auto k : speciesRangeNoI) {
            // (d / dL_0^A)
            temp = _coefs(0, 0, p + 2, k() + 1) * concentrations[_reactants[1]];
            if (volCl2 > 1) {
            for (auto j : speciesRangeNoI) {
                temp += _coefs(0, j() + 1, p + 2, k() + 1) * concentrations[_reactantMomentIds[1][j()]];
            }
            }
            values(_inverseMap(_productMomentIds[p][k()], _reactants[0]))
            += _rate(gridIndex) * temp / (double) volProd;
            // (d / dL_1^A)
            for (auto i : speciesRangeNoI) {
                temp = _coefs(i() + 1, 0, p + 2, k() + 1) * concentrations[_reactants[1]];
                if (volCl2 > 1) {
                for (auto j : speciesRangeNoI) {
                    temp += _coefs(i() + 1, j() + 1, p + 2, k() + 1) * concentrations[_reactantMomentIds[1][j()]];
                }
                }
                values(_inverseMap(_productMomentIds[p][k()], _reactantMomentIds[0][i()]))
                += _rate(gridIndex) * temp / (double) volProd;
            }
            // (d / dL_0^B)
            temp = _coefs(0, 0, p + 2, k() + 1) * concentrations[_reactants[0]];
            if (volCl1 > 1) {
            for (auto j : speciesRangeNoI) {
                temp += _coefs(j() + 1, 0, p + 2, k() + 1) * concentrations[_reactantMomentIds[0][j()]];
            }
            }
            values(_inverseMap(_productMomentIds[p][k()], _reactants[1]))
            += _rate(gridIndex) * temp / (double) volProd;
            // (d / dL_1^B)
            for (auto i : speciesRangeNoI) {
                temp = _coefs(0, i() + 1, p + 2, k() + 1) * concentrations[_reactants[0]];
                if (volCl1 > 1) {
                for (auto j : speciesRangeNoI) {
                    temp += _coefs(j() + 1, i() + 1, p + 2, k() + 1) * concentrations[_reactantMomentIds[0][j()]];
                }
                }
                values(_inverseMap(_productMomentIds[p][k()], _reactantMomentIds[1][i()]))
                += _rate(gridIndex) * temp / (double) volProd;
            }
        }
        }
    }
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::dissociationPartialDerivatives(
    ConcentrationsView concentrations,
    Kokkos::View<double*> values, std::size_t gridIndex)
{
    using AmountType = typename NetworkType::AmountType;
    constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

    // Compute the total number of elements in each cluster
    auto cl = _clusterData.getCluster(_reactants[0]);
    const auto& clReg = cl.getRegion();
    AmountType volCl = clReg.volume();
    auto prod1 = _clusterData.getCluster(_products[0]);
    const auto& prod1Reg = prod1.getRegion();
    AmountType volProd1 = prod1Reg.volume();
    auto prod2 = _clusterData.getCluster(_products[1]);
    const auto& prod2Reg = prod2.getRegion();
    AmountType volProd2 = prod2Reg.volume();

    // Compute the partials for the 0th order moments
    // First for the reactant
    double df = _rate(gridIndex) / (double) volCl;
    // Compute the values
    values(_inverseMap(_reactants[0], _reactants[0]))
    -= df * _coefs(0, 0, 0, 0);
    if (volProd1 > 1) {
        for (auto i : speciesRangeNoI) {
            values(_inverseMap(_reactants[0], _reactantMomentIds[0][i()]))
                -= df * _coefs(i() + 1, 0, 0, 0);
        }
    }
    // For the first product
    df = _rate(gridIndex) / (double) volProd1;
    values(_inverseMap(_products[0], _reactants[0])) += df * _coefs(0, 0, 0, 0);

    if (volProd1 > 1) {
        for (auto i : speciesRangeNoI) {
            values(_inverseMap(_products[0], _reactantMomentIds[0][i()]))
                += df * _coefs(i() + 1, 0, 0, 0);
        }
    }
    // For the second product
    df = _rate(gridIndex) / (double) volProd2;
    values(_inverseMap(_products[1], _reactants[0])) += df * _coefs(0, 0, 0, 0);

    if (volProd1 > 1) {
        for (auto i : speciesRangeNoI) {
            values(_inverseMap(_products[1], _reactantMomentIds[0][i()]))
                += df * _coefs(i() + 1, 0, 0, 0);
        }
    }

    // Take care of the first moments
    if (volCl > 1) {
    for (auto k : speciesRangeNoI) {
        // First for the reactant
        df = _rate(gridIndex) / (double) volCl;
        // Compute the values
        values(_inverseMap(_reactantMomentIds[0][k()], _reactants[0]))
        -= df * _coefs(0, 0, 0, k() + 1);
        for (auto i : speciesRangeNoI) {
            values(_inverseMap(_reactantMomentIds[0][k()], _reactantMomentIds[0][i()]))
            -= df * _coefs(i() + 1, 0, 0, k() + 1);
        }
        // For the first product
        if (volProd1 > 1) {
        df = _rate(gridIndex) / (double) volProd1;
        values(_inverseMap(_productMomentIds[0][k()], _reactants[0])) += df * _coefs(0, 0, 1, k() + 1);
        for (auto i : speciesRangeNoI) {
            values(_inverseMap(_productMomentIds[0][k()], _reactantMomentIds[0][i()]))
                    += df * _coefs(i() + 1, 0, 1, k() + 1);
        }
        }
        // For the second product
        if (volProd2 > 1) {
        df = _rate(gridIndex) / (double) volProd2;
        values(_inverseMap(_productMomentIds[1][k()], _reactants[0])) += df * _coefs(0, 0, 2, k() + 1);
        for (auto i : speciesRangeNoI) {
            values(_inverseMap(_productMomentIds[1][k()], _reactantMomentIds[0][i()]))
                    += df * _coefs(i() + 1, 0, 2, k() + 1);
        }
        }
    }
    }
}
}
}
