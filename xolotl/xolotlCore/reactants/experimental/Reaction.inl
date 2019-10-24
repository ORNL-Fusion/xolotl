#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TDerived>
template <typename TReactionNetwork>
Reaction<TDerived>::Reaction(TReactionNetwork& network, Type reactionType,
        std::size_t cluster0, std::size_t cluster1,
        std::size_t cluster2, std::size_t cluster3)
    :
    _type(reactionType),
    _fluxFn((_type == Type::production) ?
        &Reaction::productionFlux : &Reaction::dissociationFlux),
    _reactants((_type == Type::production) ?
        Kokkos::Array<std::size_t, 2>({cluster0, cluster1}) :
        Kokkos::Array<std::size_t, 2>({cluster0, invalid})),
    _products(_type == Type::production ?
        Kokkos::Array<std::size_t, 2>({cluster2, cluster3}) :
        Kokkos::Array<std::size_t, 2>({cluster1, cluster2})),
    _rate(static_cast<TDerived*>(this)->computeRate(network)),
    _coefs("Flux Coefficients", 5, 5, 5, 5, 5)
{
    if (_type == Type::production) {
        computeProductionCoefficients(network);
    }
    else {
        computeDissociationCoefficients(network);
    }
}


template <typename TDerived>
template <typename TReactionNetwork>
inline void
Reaction<TDerived>::computeProductionCoefficients(TReactionNetwork& network)
{
    using Species = typename TReactionNetwork::Species;
    //Compute flux coefficients
    //NOTE: _coefs is allocated 5x5x5 and zero-initialized
    //      Access with _coefs(i,j,k) with i,j,k in [0,5)
    //
    //      Get cluster tile region with, for example (plsm::Region),
    //      auto clReg = network.getCluster(_reactants[0]).getRegion()
    //
    //      Index the region with Species as aliased above (or with
    //      integers): for example
    //      auto ival = clReg[Species::V]
    //      ival represents the half-open interval (plsm::Interval)
    //      [ ival.begin(), ival.end() )
    //
    //      ival.length() is the same as ival.end() - ival.begin()
    //
    //TODO

    // Find the overlap for this reaction
    int width[4] = { };
    auto cl1Reg = network.getCluster(_reactants[0]).getRegion(),
            cl2Reg = network.getCluster(_reactants[1]).getRegion(),
            prodReg = network.getCluster(_products[0]).getRegion();
    int nOverlap = 1;
    // Loop on the phase space
    for (int i = 1; i < 5; i++) {
        // The width is the subset of the tiles for which the reaction is possible
        // For instance if we have X_1 + X_[3,5) -> X_[5,7)
        // It is only possible for 4 within X_[3,5) and 5 within X_[5,7) so the width is 1
        // More complicated with X_[3,5) + X_[5,7) -> X_[9,11)
        // 3+6, 4+5, 4+6, width is 3
        width[i - 1] = 0;
        // Loop on the interval of reactant_1 (it would be great if it could be the one with the smallest tile of the two)
        for (int j = cl1Reg[i - 1].begin(); j < cl1Reg[i - 1].end(); j++) {
            width[i - 1] += std::min(prodReg[i - 1].end()-1, cl2Reg[i - 1].end()-1 + j)
            - std::max(prodReg[i - 1].begin(), cl2Reg[i - 1].begin() + j) + 1;
        }

        nOverlap *= width[i - 1];
    }

    // Verify that overlap is larger than 0
    assert(nOverlap > 0);

    // Now we loop over the 2 dimensions of the coefs to compute all the possible sums over distances for the flux
    // The first index represents the first reactant, the second one is the second reactant
    // 0 is the 0th order, 1, 2, 3, 4, are the distances in the He, D, T, V directions (TODO make sure the order is correct)
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            if (i + j == 0) {
                // The first coefficient is simply the overlap because it is the sum over 1
                _coefs(i,j,0,0,0) = (double) nOverlap;
            }

            else if (j == 0) {
                // First order sum
                for (int l = cl2Reg[i - 1].begin(); l < cl2Reg[i - 1].end(); l++) {
                    _coefs(i,j,0,0,0) += firstOrderSum(std::max(prodReg[i - 1].begin() - l, cl1Reg[i - 1].begin()),
                            std::min(prodReg[i - 1].end()-1 - l, cl1Reg[i - 1].end()-1),
                            (double) (cl1Reg[i - 1].end()-1 + cl1Reg[i - 1].begin())/2.0);
                }
            }

            else if (i == 0) {
                // First order sum
                for (int l = cl1Reg[j - 1].begin(); l < cl1Reg[j - 1].end(); l++) {
                    _coefs(i,j,0,0,0) += firstOrderSum(std::max(prodReg[j - 1].begin() - l, cl2Reg[j - 1].begin()),
                            std::min(prodReg[j - 1].end()-1 - l, cl2Reg[j - 1].end()-1),
                            (double) (cl2Reg[j - 1].end()-1 + cl2Reg[j - 1].begin())/2.0);
                }
            }

            else {
                // Second order sum
                if (i == j) {
                    for (int l = cl1Reg[j - 1].begin(); l < cl1Reg[j - 1].end(); l++) {
                        _coefs(i,j,0,0,0) += (l - (double) (cl1Reg[j - 1].end()-1 + cl1Reg[j - 1].begin())/2.0)
                                * firstOrderSum(std::max(prodReg[j - 1].begin() - l, cl2Reg[j - 1].begin()),
                                std::min(prodReg[j - 1].end()-1 - l, cl2Reg[j - 1].end()-1),
                                (double) (cl2Reg[j - 1].end()-1 + cl2Reg[j - 1].begin())/2.0);
                    }
                }
                else {
                    _coefs(i,j,0,0,0) += _coefs(i,0,0,0,0) * _coefs(0,j,0,0,0) / (double) nOverlap;
                }
            }

            // Now we deal with the coefficients needed for the partial derivatives
            // Let's start with the product
            for (int k = 1; k < 5; k++) {

                if (i + j == 0) {
                    // First order sum
                    for (int l = cl1Reg[k - 1].begin(); l < cl1Reg[k - 1].end(); l++) {
                        _coefs(i,j,k,0,0) += firstOrderSum(std::max(prodReg[k - 1].begin(), cl2Reg[k - 1].begin() + l),
                                std::min(prodReg[k - 1].end()-1, cl2Reg[k - 1].end()-1 + l),
                                (double) (prodReg[k - 1].end()-1 + prodReg[k - 1].begin())/2.0);
                    }
                }

                else if (j == 0) {
                    // Second order sum
                    if (i == k) {
                        for (int l = cl2Reg[i - 1].begin(); l < cl2Reg[i - 1].end(); l++) {
                            _coefs(i,j,k,0,0) += secondOrderOffsetSum(std::max(prodReg[i - 1].begin() - l, cl1Reg[i - 1].begin()),
                                    std::min(prodReg[i - 1].end()-1 - l, cl1Reg[i - 1].end()-1),
                                    (double) (cl1Reg[i - 1].end()-1 + cl1Reg[i - 1].begin())/2.0,
                                    (double) (prodReg[i - 1].end()-1 + prodReg[i - 1].begin())/2.0, l);
                        }
                    }
                    else {
                        _coefs(i,j,k,0,0) += _coefs(i,0,0,0,0) * _coefs(0,0,k,0,0) / (double) nOverlap;
                    }
                }

                else if (i == 0) {
                    // Second order sum
                    if (j == k) {
                        for (int l = cl1Reg[j - 1].begin(); l < cl1Reg[j - 1].end(); l++) {
                            _coefs(i,j,k,0,0) += secondOrderOffsetSum(std::max(prodReg[j - 1].begin() - l, cl2Reg[j - 1].begin()),
                                    std::min(prodReg[j - 1].end()-1 - l, cl2Reg[j - 1].end()-1),
                                    (double) (cl2Reg[j - 1].end()-1 + cl2Reg[j - 1].begin())/2.0,
                                    (double) (prodReg[j - 1].end()-1 + prodReg[j - 1].begin())/2.0, l);
                        }
                    }
                    else {
                        _coefs(i,j,k,0,0) += _coefs(0,j,0,0,0) * _coefs(0,0,k,0,0) / (double) nOverlap;
                    }
                }

                else {
                    // Third order sum
                    if (i == j == k) {
                        for (int l = cl1Reg[i - 1].begin(); l < cl1Reg[i - 1].end(); l++) {
                            _coefs(i,j,k,0,0) += (l - (double) (cl1Reg[i - 1].end()-1 + cl1Reg[i - 1].begin())/2.0)
                                    * secondOrderOffsetSum(std::max(prodReg[i - 1].begin() - l, cl2Reg[i - 1].begin()),
                                    std::min(prodReg[i - 1].end()-1 - l, cl2Reg[i - 1].end()-1),
                                    (double) (cl2Reg[i - 1].end()-1 + cl2Reg[i - 1].begin())/2.0,
                                    (double) (prodReg[i - 1].end()-1 + prodReg[i - 1].begin())/2.0, l);
                        }
                    }
                    else if (j == k) {
                        _coefs(i,j,k,0,0) += _coefs(i,0,0,0,0) * _coefs(0,j,k,0,0) / (double) nOverlap;
                    }
                    else if (i == k) {
                        _coefs(i,j,k,0,0) += _coefs(0,j,0,0,0) * _coefs(i,0,k,0,0) / (double) nOverlap;
                    }
                    else {
                        // TODO check this is the right formula, might be divided by nOverlap^2
                        _coefs(i,j,k,0,0) += _coefs(i,0,0,0,0) * _coefs(0,j,0,0,0) * _coefs(0,0,k,0,0) / (double) nOverlap;
                    }
                }
            }
            // Let's take care of the first reactant partial derivatives
            for (int k = 1; k < 5; k++) {

                if (i + j == 0) {
                    // Same as from the flux
                    _coefs(i,j,0,k,0) += _coefs(k,0,0,0,0);
                }

                else if (j == 0) {
                    // Second order sum
                    if (i == k) {
                        for (int l = cl2Reg[i - 1].begin(); l < cl2Reg[i - 1].end(); l++) {
                            _coefs(i,j,0,k,0) += secondOrderSum(std::max(prodReg[i - 1].begin() - l, cl1Reg[i - 1].begin()),
                                    std::min(prodReg[i - 1].end()-1 - l, cl1Reg[i - 1].end()-1),
                                    (double) (cl1Reg[i - 1].end()-1 + cl1Reg[i - 1].begin())/2.0);
                        }
                    }
                    else {
                        _coefs(i,j,0,k,0) += _coefs(i,0,0,0,0) * _coefs(k,0,0,0,0) / (double) nOverlap;
                    }
                }

                else if (i == 0) {
                    _coefs(i,j,0,k,0) +=_coefs(k,j,0,0,0);
                }

                else {
                    // Third order sum
                    if (i == j == k) {
                        for (int l = cl1Reg[i - 1].begin(); l < cl1Reg[i - 1].end(); l++) {
                            _coefs(i,j,0,k,0) += (l - (double) (cl1Reg[i - 1].end()-1 + cl1Reg[i - 1].begin())/2.0)
                                    * (l - (double) (cl1Reg[i - 1].end()-1 + cl1Reg[i - 1].begin())/2.0)
                                    * firstOrderSum(std::max(prodReg[i - 1].begin() - l, cl2Reg[i - 1].begin()),
                                            std::min(prodReg[i - 1].end()-1 - l, cl2Reg[i - 1].end()-1),
                                            (double) (cl2Reg[i - 1].end()-1 + cl2Reg[i - 1].begin())/2.0);
                        }
                    }
                    else if (i == k) {
                        _coefs(i,j,0,k,0) += _coefs(0,j,0,0,0) * _coefs(i,0,0,k,0) / (double) nOverlap;
                    }
                    else if (j == k) {
                        _coefs(i,j,0,k,0) += _coefs(i,0,0,0,0) * _coefs(0,j,0,k,0) / (double) nOverlap;
                    }
                    else {
                        // TODO check this is the right formula, might be divided by nOverlap^2
                        _coefs(i,j,0,k,0) += _coefs(i,0,0,0,0) * _coefs(0,j,0,0,0) * _coefs(k,0,0,0,0) / (double) nOverlap;
                    }
                }
            }

            // Let's take care of the second reactant partial derivatives
            for (int k = 1; k < 5; k++) {

                if (i + j == 0) {
                    // Same as from the flux
                    _coefs(i,j,0,0,k) += _coefs(0,k,0,0,0);
                }

                else if (i == 0) {
                    // Second order sum
                    if (j == k) {
                        for (int l = cl1Reg[j - 1].begin(); l < cl1Reg[j - 1].end(); l++) {
                            _coefs(i,j,0,0,k) += secondOrderSum(std::max(prodReg[j - 1].begin() - l, cl2Reg[j - 1].begin()),
                                    std::min(prodReg[j - 1].end()-1 - l, cl2Reg[j - 1].end()-1),
                                    (double) (cl2Reg[j - 1].end()-1 + cl2Reg[j - 1].begin())/2.0);
                        }
                    }
                    else {
                        _coefs(i,j,0,0,k) += _coefs(0,j,0,0,0) * _coefs(0,k,0,0,0) / (double) nOverlap;
                    }
                }

                else if (j == 0) {
                    _coefs(i,j,0,0,k) +=_coefs(i,k,0,0,0);
                }

                else {
                    // Third order sum
                    if (i == j == k) {
                        for (int l = cl2Reg[i - 1].begin(); l < cl2Reg[i - 1].end(); l++) {
                            _coefs(i,j,0,0,k) += (l - (double) (cl2Reg[i - 1].end()-1 + cl2Reg[i - 1].begin())/2.0)
                                    * (l - (double) (cl2Reg[i - 1].end()-1 + cl2Reg[i - 1].begin())/2.0)
                                    * firstOrderSum(std::max(prodReg[i - 1].begin() - l, cl1Reg[i - 1].begin()),
                                            std::min(prodReg[i - 1].end()-1 - l, cl1Reg[i - 1].end()-1),
                                            (double) (cl1Reg[i - 1].end()-1 + cl1Reg[i - 1].begin())/2.0);
                        }
                    }
                    else if (i == k) {
                        _coefs(i,j,0,0,k) += _coefs(0,j,0,0,0) * _coefs(i,0,0,0,k) / (double) nOverlap;
                    }
                    else if (j == k) {
                        _coefs(i,j,0,0,k) += _coefs(i,0,0,0,0) * _coefs(0,j,0,0,k) / (double) nOverlap;
                    }
                    else {
                        // TODO check this is the right formula, might be divided by nOverlap^2
                        _coefs(i,j,0,0,k) += _coefs(i,0,0,0,0) * _coefs(0,j,0,0,0) * _coefs(0,k,0,0,0) / (double) nOverlap;
                    }
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
}
}
}
