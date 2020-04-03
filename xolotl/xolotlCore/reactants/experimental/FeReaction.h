#pragma once

#include <experimental/FeTraits.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
class FeReactionNetwork;

template <typename TSpeciesEnum>
class FeProductionReaction :
    public ProductionReaction<FeReactionNetwork<TSpeciesEnum>,
        FeProductionReaction<TSpeciesEnum>>
{
public:
    using Superclass = ProductionReaction<FeReactionNetwork<TSpeciesEnum>,
        FeProductionReaction<TSpeciesEnum>>;

    using Superclass::Superclass;
};

template <typename TSpeciesEnum>
class FeDissociationReaction :
    public DissociationReaction<FeReactionNetwork<TSpeciesEnum>,
        FeDissociationReaction<TSpeciesEnum>>
{
public:
    using Superclass = DissociationReaction<FeReactionNetwork<TSpeciesEnum>,
            FeDissociationReaction<TSpeciesEnum>>;

    using Superclass::Superclass;

    KOKKOS_INLINE_FUNCTION
    double
    computeBindingEnergy();
};

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
FeDissociationReaction<TSpeciesEnum>::computeBindingEnergy()
{
    using NetworkType = typename Superclass::NetworkType;
    using Species = typename NetworkType::Species;
    using Composition = typename NetworkType::Composition;

    constexpr double heTrapTable[9] = {
        0.0, 4.31, 2.90, 2.02, 1.09, 0.58, 0.13, -0.25, -0.59
    };

    double be = 5.0;

    auto cl = this->_clusterData.getCluster(this->_reactant);
    auto prod1 = this->_clusterData.getCluster(this->_products[0]);
    auto prod2 = this->_clusterData.getCluster(this->_products[1]);

    auto clReg = cl.getRegion();
    auto prod1Reg = prod1.getRegion();
    auto prod2Reg = prod2.getRegion();
    bool useTable = false;
    if (clReg.isSimplex() && prod1Reg.isSimplex() && prod2Reg.isSimplex()) {
        Composition comp = clReg.getOrigin();
        Composition prod1Comp = prod1Reg.getOrigin();
        Composition prod2Comp = prod2Reg.getOrigin();
        if (comp.isOnAxis(Species::He)) {
            if (prod1Comp.isOnAxis(Species::He) || prod2Comp.isOnAxis(Species::He)) {
                if (comp[Species::He] == 2)
                    be = 0.5;
                else
                    be = 1.0;
            }
            if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
                be = heTrapTable[comp[Species::He]];
            }
        }
        else if (comp.isOnAxis(Species::V)) {
            auto size = comp[Species::V];
            be = 1.73 - 2.59
                * (pow((double) size, 2.0 / 3.0)
                - pow((double) size - 1.0, 2.0 / 3.0));
        }
        else if (comp.isOnAxis(Species::I)) {
            // Nothing
        }
        else {
            // HeV
            if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
                auto amtHe = comp[Species::He], amtV = comp[Species::V];
                be = 1.73 - 2.59
                    * (pow((double) amtV, 2.0 / 3.0)
                    - pow((double) amtV - 1.0, 2.0 / 3.0))
                    + 2.5 * log(1.0 + ((double) amtHe / (double) amtV));
            }
            if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
                auto amtHe = comp[Species::He], amtV = comp[Species::V];
                be = 4.88 + 2.59
                    * (pow((double) amtV, 2.0 / 3.0)
                    - pow((double) amtV - 1.0, 2.0 / 3.0))
                    - 2.5 * log(1.0 + ((double) amtHe / (double) amtV));
            }
        }
    }
    else {
        Composition lo = clReg.getOrigin();
        Composition hi = clReg.getUpperLimitPoint();
        Composition prod1Comp = prod1Reg.getOrigin();
        Composition prod2Comp = prod2Reg.getOrigin();
        // HeV
        if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
            double amtHe = (double) (lo[Species::He] + hi[Species::He] - 1) / 2.0,
                amtV = (double) (lo[Species::V] + hi[Species::V] - 1) / 2.0;
            be = 1.73 - 2.59
                * (pow((double) amtV, 2.0 / 3.0)
                - pow((double) amtV - 1.0, 2.0 / 3.0))
                + 2.5 * log(1.0 + ((double) amtHe / (double) amtV));
        }
        if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
            double amtHe = (double) (lo[Species::He] + hi[Species::He] - 1) / 2.0,
                amtV = (double) (lo[Species::V] + hi[Species::V] - 1) / 2.0;
            be = 4.88 + 2.59
                * (pow((double) amtV, 2.0 / 3.0)
                - pow((double) amtV - 1.0, 2.0 / 3.0))
                - 2.5 * log(1.0 + ((double) amtHe / (double) amtV));
        }
    }

    return max(be, -5.0);
}
}
}
