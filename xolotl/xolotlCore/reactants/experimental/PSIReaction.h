#pragma once

#include <experimental/PSITraits.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
class PSIReactionNetwork;

template <typename TSpeciesEnum>
class PSIProductionReaction :
    public ProductionReaction<PSIReactionNetwork<TSpeciesEnum>,
        PSIProductionReaction<TSpeciesEnum>>
{
public:
    using Superclass = ProductionReaction<PSIReactionNetwork<TSpeciesEnum>,
        PSIProductionReaction<TSpeciesEnum>>;

    using Superclass::Superclass;
};

template <typename TSpeciesEnum>
class PSIDissociationReaction :
    public DissociationReaction<PSIReactionNetwork<TSpeciesEnum>,
        PSIDissociationReaction<TSpeciesEnum>>
{
public:
    using Superclass = DissociationReaction<PSIReactionNetwork<TSpeciesEnum>,
            PSIDissociationReaction<TSpeciesEnum>>;

    using Superclass::Superclass;

    KOKKOS_INLINE_FUNCTION
    double
    computeBindingEnergy();
};

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSIDissociationReaction<TSpeciesEnum>::computeBindingEnergy()
{
    using NetworkType = typename Superclass::NetworkType;

    constexpr double beTableV1[10][7] = {
        //H:  1     2     3     4     5     6      // He:
        {0.0, 1.21, 1.17, 1.05, 0.93, 0.85, 0.60}, // 0
        {0.0, 1.00, 0.95, 0.90, 0.88, 0.80, 0.60}, // 1
        {0.0, 0.96, 0.92, 0.85, 0.84, 0.83, 0.50}, // 2
        {0.0, 0.86, 0.81, 0.69, 0.64, 0.65, 0.50}, // 3
        {0.0, 0.83, 0.80, 0.65, 0.60, 0.60, 0.55}, // 4
        {0.0, 0.83, 0.80, 0.60, 0.50, 0.50, 0.50}, // 5
        {0.0, 0.80, 0.70, 0.60, 0.50, 0.50, 0.50}, // 6
        {0.0, 0.80, 0.75, 0.65, 0.55, 0.55, 0.45}, // 7
        {0.0, 0.80, 0.80, 0.70, 0.65, 0.60, 0.55}, // 8
        {0.0, 0.80, 0.80, 0.75, 0.70, 0.65, 0.60}, // 9
    };

    static constexpr double beTableV2[15][12] = {
        //H:  1     2     3     4     5     6     7     8     9     10    11     // He:
        {0.0, 1.63, 1.31, 1.25, 1.16, 1.00, 1.00, 0.95, 0.95, 0.75, 0.70, 0.65}, // 0
        {0.0, 1.30, 1.30, 1.24, 1.08, 0.95, 0.95, 0.95, 0.95, 0.75, 0.70, 0.65}, // 1
        {0.0, 1.15, 1.14, 1.11, 1.14, 0.95, 0.95, 0.95, 0.90, 0.75, 0.70, 0.65}, // 2
        {0.0, 1.12, 1.06, 0.99, 0.99, 0.90, 0.95, 0.90, 0.90, 0.70, 0.70, 0.65}, // 3
        {0.0, 1.10, 1.06, 0.99, 0.99, 0.90, 0.95, 0.90, 0.90, 0.70, 0.65, 0.65}, // 4
        {0.0, 1.10, 1.05, 0.99, 0.99, 0.90, 0.90, 0.90, 0.90, 0.70, 0.65, 0.65}, // 5
        {0.0, 1.10, 1.05, 0.99, 0.99, 0.90, 0.90, 0.90, 0.85, 0.70, 0.65, 0.60}, // 6
        {0.0, 1.05, 1.00, 0.95, 0.95, 0.90, 0.90, 0.90, 0.85, 0.65, 0.65, 0.60}, // 7
        {0.0, 1.05, 1.00, 0.95, 0.95, 0.90, 0.90, 0.85, 0.85, 0.65, 0.65, 0.60}, // 8
        {0.0, 1.05, 1.00, 0.95, 0.95, 0.85, 0.85, 0.85, 0.85, 0.65, 0.65, 0.60}, // 9
        {0.0, 1.00, 0.95, 0.90, 0.90, 0.85, 0.85, 0.85, 0.80, 0.65, 0.60, 0.60}, // 10
        {0.0, 0.95, 0.95, 0.90, 0.90, 0.85, 0.85, 0.85, 0.80, 0.65, 0.60, 0.60}, // 11
        {0.0, 0.95, 0.90, 0.90, 0.85, 0.85, 0.85, 0.80, 0.80, 0.60, 0.60, 0.55}, // 12
        {0.0, 0.90, 0.90, 0.85, 0.85, 0.85, 0.85, 0.80, 0.80, 0.60, 0.60, 0.55}, // 13
        {0.0, 0.90, 0.90, 0.85, 0.85, 0.80, 0.80, 0.80, 0.70, 0.60, 0.60, 0.55}, // 14
    };

    using Species = typename NetworkType::Species;
    using Composition = typename NetworkType::Composition;

    double be = 0.0;

    auto cl = this->_clusterData.getCluster(this->_reactant);
    auto prod1 = this->_clusterData.getCluster(this->_products[0]);
    auto prod2 = this->_clusterData.getCluster(this->_products[1]);

    auto clReg = cl.getRegion();
    auto prod1Reg = prod1.getRegion();
    auto prod2Reg = prod2.getRegion();
    bool useTable = false;
    if (clReg.isSimplex()) {
        if (prod1Reg.isSimplex()) {
            auto orig1 = prod1Reg.getOrigin();
            if (orig1.isOnAxis(Species::D) || orig1.isOnAxis(Species::T)) {
                useTable = true;
            }
        }
        else if (prod2Reg.isSimplex()) {
            auto orig2 = prod2Reg.getOrigin();
            if (orig2.isOnAxis(Species::D) || orig2.isOnAxis(Species::T)) {
                useTable = true;
            }
        }
    }

    if (useTable) {
        Composition comp(clReg.getOrigin());
        auto hAmount = comp[Species::D] + comp[Species::T];
        if (comp[Species::V] == 1) {
            be = beTableV1[comp[Species::He]][hAmount];
        }
        else if (comp[Species::V] == 2) {
            be = beTableV2[comp[Species::He]][hAmount];
        }
    }

    if (be == 0.0) {
        be = prod1.getFormationEnergy() + prod2.getFormationEnergy() -
            cl.getFormationEnergy();
    }

    return max(be, -5.0);
}
}
}
