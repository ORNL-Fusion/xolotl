#pragma once

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
NEReactionGenerator::operator()(IndexType i, IndexType j, TTag tag) const
{
    using Species = typename Network::Species;
    using Composition = typename Network::Composition;
    using AmountType = typename Network::AmountType;

    auto numClusters = this->getNumberOfClusters();

    // Get the composition of each cluster
    const auto& cl1Reg = this->getCluster(i).getRegion();
    const auto& cl2Reg = this->getCluster(j).getRegion();
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
    for (IndexType k = 0; k < numClusters; ++k) {
        // Get the composition
        const auto& prodReg = this->getCluster(k).getRegion();
        // Check the bounds
        if (prodReg[Species::Xe].begin() > bounds.second) {
            continue;
        }
        else if (prodReg[Species::Xe].end() - 1 < bounds.first) {
            continue;
        }

        this->addProductionReaction(tag, {i, j, k});

        if (!cl1Reg.isSimplex() && !cl2Reg.isSimplex()) {
            continue;
        }
        // Is the size of one of them one?
        if (lo1[Species::Xe] == 1 || lo2[Species::Xe] == 1) {
            this->addDissociationReaction(tag, {k, i, j});
        }
    }
}

template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
NEReactionGenerator::addSinks(IndexType i, TTag tag) const
{
    // Nothing
}
}

inline
detail::NEReactionGenerator
NEReactionNetwork::getReactionGenerator() const noexcept
{
    return detail::NEReactionGenerator{*this};
}
}
}
