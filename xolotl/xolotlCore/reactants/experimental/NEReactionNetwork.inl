#pragma once

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
KOKKOS_INLINE_FUNCTION
void
NEReactionValidator::operator()(std::size_t i, std::size_t j,
    const Subpaving& subpaving, const UpperTriangle<ClusterSet>& prodSet,
    const UpperTriangle<ClusterSet>& dissSet) const
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

    // General case
    Kokkos::pair<AmountType, AmountType> bounds;
    // Compute the bounds
    auto low = lo1[Species::Xe] + lo2[Species::Xe];
    auto high = hi1[Species::Xe] + hi2[Species::Xe] - 2;
    bounds = {low, high};

    // Look for potential product
    for (std::size_t k = 0; k < numClusters; ++k) {
        // Get the composition
        const auto& prodReg = tiles(k).getRegion();
        bool isGood = true;
        // Check the bounds
        if (prodReg[Species::Xe].begin() > bounds.second) {
            isGood = false;
        }
        else if (prodReg[Species::Xe].end() - 1 < bounds.first) {
            isGood = false;
        }

        if (isGood) {
            // Increase nProd
            prodSet(i, j) = {i, j, k};
            // TODO: will have to add some rules, i or j should be a simplex cluster of max size 1
            if (!cl1Reg.isSimplex() && !cl2Reg.isSimplex()) {
                continue;
            }
            // Is the size of one of them one?
            if (lo1[Species::Xe] == 1 || lo2[Species::Xe] == 1) {
                dissSet(i, j) = {k, i, j};
            }
        }
    }
}
}

inline
detail::NEReactionValidator
NEReactionNetwork::getReactionValidator() const noexcept
{
    return {};
}
}
}
