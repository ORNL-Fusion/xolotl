#pragma once

#include <xolotl/core/network/AlloyTraits.h>
#include <xolotl/core/network/SinkReaction.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
class AlloyReactionNetwork;

class AlloyProductionReaction :
    public ProductionReaction<AlloyReactionNetwork, AlloyProductionReaction>
{
public:
    using Superclass =
        ProductionReaction<AlloyReactionNetwork, AlloyProductionReaction>;

    using Superclass::Superclass;
};

class AlloyDissociationReaction :
    public DissociationReaction<AlloyReactionNetwork, AlloyDissociationReaction>
{
public:
    using Superclass =
        DissociationReaction<AlloyReactionNetwork, AlloyDissociationReaction>;

    using Superclass::Superclass;

    KOKKOS_INLINE_FUNCTION
    double
    computeBindingEnergy()
    {
        using Species = typename Superclass::Species;
        using Composition = typename Superclass::Composition;

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
            if (comp.isOnAxis(Species::Void)) {
            	double n = comp[Species::Void];
                if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
                	be = 3.5 - 3.45 * (pow(n + 1.0, 2.0 / 3.0) - pow(n, 2.0 / 3.0));
                }
                else if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
                    be = 1.9 - 3.1 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
                }
            }
            else if (comp.isOnAxis(Species::Faulted)) {
            	double n = comp[Species::Faulted];
                if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
                	be = 1.9 - 3.2 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
                }
            }
            else if (comp.isOnAxis(Species::V)) {
            	double n = comp[Species::V];
                if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
                	be = 1.9 - 3.1 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
                }
            }
            else if (comp.isOnAxis(Species::I)) {
            	double n = comp[Species::I];
                if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
                	be = 3.5 - 2.5 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
                }
            }
        }
//        else {
//            Composition lo = clReg.getOrigin();
//            Composition hi = clReg.getUpperLimitPoint();
//            Composition prod1Comp = prod1Reg.getOrigin();
//            Composition prod2Comp = prod2Reg.getOrigin();
//            // HeV
//            double amtHe = (double) (lo[Species::He] + hi[Species::He] - 1) / 2.0,
//                amtV = (double) (lo[Species::V] + hi[Species::V] - 1) / 2.0;
//            if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
//                be = 1.73 - 2.59
//                    * (pow(amtV, 2.0 / 3.0)
//                    - pow(amtV - 1.0, 2.0 / 3.0))
//                    + 2.5 * log(1.0 + (amtHe / amtV));
//            }
//            if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
//                be = 4.88 + 2.59
//                    * (pow(amtV, 2.0 / 3.0)
//                    - pow(amtV - 1.0, 2.0 / 3.0))
//                    - 2.5 * log(1.0 + (amtHe / amtV));
//            }
//        }

        return util::max(0.1, be);
    }
};

class AlloySinkReaction :
    public SinkReaction<AlloyReactionNetwork, AlloySinkReaction>
{
public:
    using Superclass =
        SinkReaction<AlloyReactionNetwork, AlloySinkReaction>;

    using Superclass::Superclass;

    KOKKOS_INLINE_FUNCTION
    double
    getSinkBias()
    {
        using Species = typename Superclass::Species;
        using Composition = typename Superclass::Composition;

        double bias = 1.0;

        auto cl = this->_clusterData.getCluster(this->_reactant);

        auto clReg = cl.getRegion();
        if (clReg.isSimplex()) {
            Composition comp = clReg.getOrigin();
            if (comp.isOnAxis(Species::I)) {
                bias = 1.2;
            }
        }

        return bias;
    }
};
}
}
}
