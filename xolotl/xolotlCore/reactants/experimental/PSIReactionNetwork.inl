#pragma once

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
void
PSIReactionGenerator<TSpeciesEnum>::operator()(std::size_t i, std::size_t j,
    const Subpaving& subpaving, const UpperTriangle<Kokkos::pair<ClusterSet, ClusterSet> >& prodSet,
    const UpperTriangle<Kokkos::pair<ClusterSet, ClusterSet> >& dissSet) const
{
    using Species = typename Network::Species;
    using Composition = typename Network::Composition;
    using AmountType = typename Network::AmountType;

    constexpr auto species = Network::getSpeciesRange();
    constexpr auto speciesNoI = Network::getSpeciesRangeNoI();
    constexpr auto invalid = Network::invalid;

    const auto& tiles = subpaving.getTiles(plsm::onDevice);
    auto numClusters = tiles.extent(0);

    // Get the composition of each cluster
    const auto& cl1Reg = tiles(i).getRegion();
    const auto& cl2Reg = tiles(j).getRegion();
    Composition lo1 = cl1Reg.getOrigin();
    Composition hi1 = cl1Reg.getUpperLimitPoint();
    Composition lo2 = cl2Reg.getOrigin();
    Composition hi2 = cl2Reg.getUpperLimitPoint();

    // Special case for I + I
    if (cl1Reg.isSimplex() && cl2Reg.isSimplex() && lo1.isOnAxis(Species::I)
        && lo2.isOnAxis(Species::I)) {
        // Compute the composition of the new cluster
        auto size = lo1[Species::I] + lo2[Species::I];
        // Find the corresponding cluster
        Composition comp = Composition::zero();
        comp[Species::I] = size;
        auto iProdId = subpaving.findTileId(comp, plsm::onDevice);
        if (iProdId != invalid) {
            prodSet(i, j).first = {i, j, iProdId};
            if (lo1[Species::I] == 1 || lo2[Species::I] == 1) {
                dissSet(i, j).first = {iProdId, i, j};
            }
        }
        return;
    }

    // Special case for I + V
    if (cl1Reg.isSimplex() && cl2Reg.isSimplex() &&
            ((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
             (lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::I)))) {
        // Find out which one is which
        auto vSize = lo1.isOnAxis(Species::V)? lo1[Species::V] :
           lo2[Species::V];
        auto iSize = lo1.isOnAxis(Species::I)? lo1[Species::I] :
           lo2[Species::I];
        // Compute the product size
        int prodSize = vSize - iSize;
        // 3 cases
        if (prodSize > 0) {
            // Looking for V cluster
            Composition comp = Composition::zero();
            comp[Species::V] = prodSize;
            auto vProdId = subpaving.findTileId(comp, plsm::onDevice);
            if (vProdId != invalid) {
                prodSet(i, j).first = {i, j, vProdId};
                // No dissociation
            }
        }
        else if (prodSize < 0) {
            // Looking for I cluster
            Composition comp = Composition::zero();
            comp[Species::I] = -prodSize;
            auto iProdId = subpaving.findTileId(comp, plsm::onDevice);
            if (iProdId != invalid) {
                prodSet(i, j).first = {i, j, iProdId};
                // No dissociation
            }
        }
        else {
            // No product
            prodSet(i, j).first = {i, j};
        }
        return;
    }

    // General case
    constexpr auto numSpeciesNoI = Network::getNumberOfSpeciesNoI();
    using BoundsArray =
        Kokkos::Array<Kokkos::pair<AmountType, AmountType>, numSpeciesNoI>;
    plsm::EnumIndexed<BoundsArray, Species> bounds;
    // Loop on the species
    for (auto l : species) {
        auto low = lo1[l] + lo2[l];
        auto high = hi1[l] + hi2[l] - 2;
        // Special case for I
        if (l == Species::I) {
            bounds[Species::V].first -= high;
            bounds[Species::V].second -= low;
        }
        else {
            bounds[l] = {low, high};
        }
    }

    // Look for potential product
    std::size_t nProd = 0;
    for (std::size_t k = 0; k < numClusters; ++k) {
        // Get the composition
        const auto& prodReg = tiles(k).getRegion();
        bool isGood = true;
        // Loop on the species
        // TODO: check l correspond to the same species in bounds and prod
        for (auto l : speciesNoI) {
            if (prodReg[l()].begin() > bounds[l()].second) {
                isGood = false;
                break;
            }
            if (prodReg[l()].end() - 1 < bounds[l()].first) {
                isGood = false;
                break;
            }
        }

        if (isGood) {
            // Increase nProd
            nProd++;
            prodSet(i, j).first = {i, j, k};
            // TODO: will have to add some rules, i or j should be a simplex cluster of max size 1
            if (!cl1Reg.isSimplex() && !cl2Reg.isSimplex()) {
                continue;
            }
            // Loop on the species
            bool isOnAxis1 = false, isOnAxis2 = false;
            for (auto l : species) {
                if (lo1.isOnAxis(l()) && lo1[l()] == 1) isOnAxis1 = true;
                if (lo2.isOnAxis(l()) && lo2[l()] == 1) isOnAxis2 = true;
            }
            if (isOnAxis1 || isOnAxis2) {
                if (lo1.isOnAxis(Species::I) || lo2.isOnAxis(Species::I)) {
                    continue;
                }

                dissSet(i, j).first = {k, i, j};
            }
        }
    }

    // Special case for trap-mutation
    if (nProd == 0) {
        // Look for larger clusters only if one of the reactant is pure He
        if (!(cl1Reg.isSimplex() && lo1.isOnAxis(Species::He)) &&
                !(cl2Reg.isSimplex() && lo2.isOnAxis(Species::He))) {
            return;
        }

        // Check that both reactants contain He
        if (cl1Reg[Species::He].begin() < 1 ||
                cl2Reg[Species::He].begin() < 1) {
            return;
        }

        // Loop on possible I sizes
        // TODO: get the correct value for maxISize
        std::size_t maxISize = 6;
        for (std::size_t n = 1; n <= maxISize; ++n) {
            // Find the corresponding cluster
            Composition comp = Composition::zero();
            comp[Species::I] = n;
            auto iClusterId = subpaving.findTileId(comp, plsm::onDevice);

            bounds[Species::V].first += 1;
            bounds[Species::V].second += 1;

            // Look for potential product
            std::size_t nProd = 0;
            for (std::size_t k = 0; k < numClusters; ++k) {
                // Get the composition
                const auto& prodReg = tiles(k).getRegion();
                bool isGood = true;
                // Loop on the species
                // TODO: check l correspond to the same species in bounds and prod
                for (auto l : speciesNoI) {
                    if (prodReg[l()].begin() > bounds[l()].second) {
                        isGood = false;
                        break;
                    }
                    if (prodReg[l()].end() - 1 < bounds[l()].first) {
                        isGood = false;
                        break;
                    }
                }

                if (isGood) {
                    // Increase nProd
                    nProd++;
                    prodSet(i, j).first = {i, j, k, iClusterId};
                    // No dissociation
                }
            }
            // Stop if we found a product
            if (nProd > 0) {
                break;
            }
        }
    }
}
}
}
}
