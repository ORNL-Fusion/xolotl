#pragma once

#include <experimental/ReactionNetwork.h>
#include <experimental/FeReaction.h>
#include <experimental/FeTraits.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
class FeReactionGenerator;
}

class FeReactionNetwork :
    public ReactionNetwork<FeReactionNetwork>
{
    friend class ReactionNetwork<FeReactionNetwork>;
    friend class detail::ReactionNetworkWorker<FeReactionNetwork>;

public:
    using Superclass = ReactionNetwork<FeReactionNetwork>;
    using Subpaving = typename Superclass::Subpaving;
    using Composition = typename Superclass::Composition;
    using Species = typename Superclass::Species;
    using AmountType = typename Superclass::AmountType;
    using IndexType = typename Superclass::IndexType;
    using ConcentrationsView = typename Superclass::ConcentrationsView;
    using FluxesView = typename Superclass::FluxesView;
    using ClusterDataRef = typename Superclass::ClusterDataRef;

    using Superclass::Superclass;

    void
    checkTiles(const IOptions& options)
    {
        auto maxHe = static_cast<AmountType>(options.getMaxImpurity() + 1);
        auto maxV = static_cast<AmountType>(options.getMaxV() + 1);
        auto& tiles = this->getSubpaving().getTiles(plsm::onDevice);
        auto numClusters = tiles.extent(0);
        Kokkos::parallel_for(numClusters, KOKKOS_LAMBDA (const IndexType i) {
            auto clReg = tiles(i).getRegion();
            if (clReg[Species::He].end() > maxHe) {
                Region r({Ival{clReg[Species::He].begin(), maxHe},
                    Ival{clReg[Species::V].begin(), clReg[Species::V].end()},
                    Ival{clReg[Species::I].begin(), clReg[Species::I].end()}});
                auto id = tiles(i).getOwningZoneIndex();
                auto newTile = plsm::Tile<Region>(r, id);
                tiles(i) = newTile;
            }
            clReg = tiles(i).getRegion();
            if (clReg[Species::He].begin() == 0 && clReg[Species::V].begin() > 10) {
                Region r({Ival{1, clReg[Species::He].end()},
                    Ival{clReg[Species::V].begin(), clReg[Species::V].end()},
                    Ival{clReg[Species::I].begin(), clReg[Species::I].end()}});
                auto id = tiles(i).getOwningZoneIndex();
                auto newTile = plsm::Tile<Region>(r, id);
                tiles(i) = newTile;
            }
            clReg = tiles(i).getRegion();
            if (clReg[Species::V].end() > maxV) {
                Region r({Ival{clReg[Species::He].begin(), clReg[Species::He].end()},
                    Ival{clReg[Species::V].begin(), maxV},
                    Ival{clReg[Species::I].begin(), clReg[Species::I].end()}});
                auto id = tiles(i).getOwningZoneIndex();
                auto newTile = plsm::Tile<Region>(r, id);
                tiles(i) = newTile;
            }
            clReg = tiles(i).getRegion();
            if (clReg[Species::V].begin() == 0 && clReg[Species::He].begin() > 8) {
                Region r({Ival{clReg[Species::He].begin(), clReg[Species::He].end()},
                    Ival{1, clReg[Species::V].end()},
                    Ival{clReg[Species::I].begin(), clReg[Species::I].end()}});
                auto id = tiles(i).getOwningZoneIndex();
                auto newTile = plsm::Tile<Region>(r, id);
                tiles(i) = newTile;
            }
        });
    }

private:
    double
    checkLatticeParameter(double latticeParameter)
    {
        if (latticeParameter <= 0.0) {
            return ironLatticeConstant;
        }
        return latticeParameter;
    }

    double
    checkImpurityRadius(double impurityRadius)
    {
        if (impurityRadius <= 0.0) {
            return heliumRadius;
        }
        return impurityRadius;
    }

    detail::FeReactionGenerator
    getReactionGenerator() const noexcept;
};

namespace detail
{
class FeReactionGenerator :
    public ReactionGenerator<FeReactionNetwork, FeReactionGenerator>
{
    friend class ReactionGeneratorBase<FeReactionNetwork, FeReactionGenerator>;

public:
    using NetworkType = FeReactionNetwork;
    using Subpaving = typename NetworkType::Subpaving;
    using IndexType = typename NetworkType::IndexType;

    using Superclass = ReactionGenerator<FeReactionNetwork,
        FeReactionGenerator>;

    using Superclass::Superclass;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    operator()(IndexType i, IndexType j, TTag tag) const;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    addSinks(IndexType i, TTag tag) const;

private:
    ReactionCollection<NetworkType>
    getReactionCollection() const;
};
}
}
}

#include <experimental/FeClusterGenerator.h>
#include <experimental/FeReactionNetwork.inl>
