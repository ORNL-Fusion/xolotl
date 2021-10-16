#pragma once

#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TRegion>
KOKKOS_INLINE_FUNCTION
double
getRate(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg, const double r0,
	const double r1, const double dc0, const double dc1)
{
	//constexpr double pi = ::xolotl::core::pi;
	//constexpr double rCore = ::xolotl::core::alphaZrCoreRadius;
	//const double zs = 4.0 * pi * (r0 + r1 + rCore);

	//using Species = typename TRegion::EnumIndex;
	//xolotl::core::network::detail::Composition<typename TRegion::VectorType,
	//	Species>
	//	lo0 = pairCl0Reg.getOrigin();
	//xolotl::core::network::detail::Composition<typename TRegion::VectorType,
	//	Species>
	//	lo1 = pairCl1Reg.getOrigin();
	//bool cl0Is1D = (lo0[Species::I] == 9), cl1Is1D = (lo1[Species::I] == 9);

    // Until the macroscopic cross-section method is established, just assume 3-D rates of reactions:
	// Cluster 0 is 1D diffuser
	//if (cl0Is1D) {
	//	return zs * (dc0 + dc1);
	//}

	// Cluster 1 is 1D diffuser
	//else if (cl1Is1D) {
	//	return zs * (dc0 + dc1);
	//}

	// None of them is a 1D diffuser
	//return zs * (dc0 + dc1);


    constexpr double pi = ::xolotl::core::pi;
    constexpr double rCore = ::xolotl::core::alphaZrCoreRadius;
    const double zs = 4.0 * pi * (r0 + r1 + rCore);

    using Species = typename TRegion::EnumIndex;
    xolotl::core::network::detail::Composition<typename TRegion::VectorType,
        Species>
        lo0 = pairCl0Reg.getOrigin();
    xolotl::core::network::detail::Composition<typename TRegion::VectorType,
        Species>
        lo1 = pairCl1Reg.getOrigin();
    bool cl0Is1D = (lo0[Species::I] == 9), cl1Is1D = (lo1[Species::I] == 9);

    //Determine if clusters are vacancy or interstitial and initialize variables
    bool cl0IsV = lo0.isOnAxis(Species::V);
    bool cl1IsV = lo1.isOnAxis(Species::V);
    double n0 = 0; //size of cluster 0
    double n1 = 0; //size of cluster 1
    double rdI = 6.0; //dislocation capture radius for I defects
    double rdV = 3.0; //dislocation capture radius for V defects
    double rd = 0.0; //dislocation capture radius
    double p = 1.0; //anisotropy ratio for diffusing defect
    double Pl = 1.0; //Capture efficiency for diffusing defect

    //Determine parameters based on cluster type and size
    if (cl0IsV) n0 = lo0[Species::V];
    else n0 = lo0[Species::I];
    bool cl0IsLoop = (n0 > 9);

    if (cl1IsV) double n1 = lo1[Species::V];
    else double n1 = lo1[Species::I];
    bool cl1IsLoop = (n1 > 9);

    //Need to ask Sophie how to access anisotropy ratio, and how to add capture efficiencies to extraData
    //For now, manually assign anisotropy ratios for 573 K
    if ((n0 == 1 && cl0IsV) || (n1 == 1 && cl1IsV)) p = 0.80; //Mobile vacancy, n = 1
    else if ((n0 == 1 && !cl0IsV) || (n1 == 1 && !cl1IsV)) p = 0.77; //Mobile interstitial, n = 1
    if ((n0 == 2 && cl0IsV) || (n1 == 2 && cl1IsV)) p = 1.17;
    else if ((n0 == 2 && !cl0IsV) || (n1 == 2 && !cl1IsV)) p = 0.50;
    if ((n0 == 3 && cl0IsV) || (n1 == 3 && cl1IsV)) p = 0.93;
    else if ((n0 == 3 && !cl0IsV) || (n1 == 3 && !cl1IsV)) p = 0.35;
    if ((n0 == 4 && cl0IsV) || (n1 == 4 && cl1IsV)) p = 0.39;
    else if ((n0 == 4 && !cl0IsV) || (n1 == 4 && !cl1IsV)) p = 0.31;
    if ((n0 == 5 && cl0IsV) || (n1 == 5 && cl1IsV)) p = 3.17;
    else if ((n0 == 5 && !cl0IsV) || (n1 == 5 && !cl1IsV)) p = 0.19;
    if ((n0 == 6 && cl0IsV) || (n1 == 6 && cl1IsV)) p = 1.02;
    if ((n0 == 7 && cl0IsV) || (n1 == 7 && cl1IsV)) p = 1.45;
    if ((n0 == 8 && cl0IsV) || (n1 == 8 && cl1IsV)) p = 1.49;
    if ((n0 == 9 && cl0IsV) || (n1 == 9 && cl1IsV)) p = 0.84;
    else if ((n0 == 9 && !cl0IsV) || (n1 == 9 && !cl1IsV)) p = 0.50;

    // These rates are only for 3-D mobile diffusers, but for now assume that it works for 1-D diffusers as well
    // Cluster 0 is a dislocation loop:
    if (cl0IsLoop) {
        if (cl1IsV) rd = rdV;
        else rd = rdI;
        double alpha = pow(1+pow(r0 / (3*(r1 + rd)), 2), -1);
        double rateSpherical = 4.0 * pi * (r0 + r1 + rd);
        double rateToroidal = (4.0 * pi * pi * r0) / log(1 + (8 * r0) / (r1 + rd));

        //For now, manually assign the capture efficiency (assuming only prismatic clusters)
        if (cl0IsV) Pl = 0.78 * pow(p,-2) + 0.66 * p - 0.44;
        else Pl = 0.70 * pow(p,-2) + 0.78 * p - 0.47;

        return (1 - alpha) * rateToroidal * Pl + alpha * rateSpherical;
    }

    // Cluster 1 is a dislocation loop:
    if (cl1IsLoop) {
        if (cl0IsV) rd = rdV;
        else rd = rdI;
        double alpha = pow(1+pow(r1 / (3*(r0 + rd)), 2), -1);
        double rateSpherical = 4.0 * pi * (r0 + r1 + rd);
        double rateToroidal = (4.0 * pi * pi * r1) / log(1 + (8 * r1) / (r0 + rd));

        //For now, manually assign the capture efficiency (assuming only prismatic clusters)
        if (cl1IsV) Pl = 0.78 * pow(p,-2) + 0.66 * p - 0.44;
        else Pl = 0.70 * pow(p,-2) + 0.78 * p - 0.47;

        return (1 - alpha) * rateToroidal * Pl + alpha * rateSpherical;
    }

    // None of the clusters are loops (interaction is based on spherical volume)
    return zs * (dc0 + dc1);
    }

KOKKOS_INLINE_FUNCTION
double
ZrProductionReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

KOKKOS_INLINE_FUNCTION
double
ZrDissociationReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

KOKKOS_INLINE_FUNCTION
double
ZrDissociationReaction::computeBindingEnergy()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	double be = 5.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);
	auto prod1 = this->_clusterData->getCluster(this->_products[0]);
	auto prod2 = this->_clusterData->getCluster(this->_products[1]);

	auto clReg = cl.getRegion();
	auto prod1Reg = prod1.getRegion();
	auto prod2Reg = prod2.getRegion();
	Composition lo = clReg.getOrigin();
	Composition hi = clReg.getUpperLimitPoint();
	Composition prod1Comp = prod1Reg.getOrigin();
	Composition prod2Comp = prod2Reg.getOrigin();

	if (lo.isOnAxis(Species::V)) {
		double n = (double)(lo[Species::V] + hi[Species::V] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			if (n < 18) be = 2.03 - 1.9 * (pow(n, 0.84) - pow(n - 1.0, 0.84));
            else be = 2.03 - 3.4 * (pow(n, 0.70) - pow(n - 1.0, 0.70));

		}
	}
	else if (lo.isOnAxis(Species::I)) {
		double n = (double)(lo[Species::I] + hi[Species::I] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
			if (n < 7) be = 2.94 - 2.8 * (pow(n, 0.81) - pow(n - 1.0, 0.81));
            else be = 2.94 - 4.6 * (pow(n, 0.66) - pow(n - 1.0, 0.66));

		}
	}

	return util::max(0.1, be);
}

KOKKOS_INLINE_FUNCTION
double
ZrSinkReaction::computeRate(IndexType gridIndex)
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	auto cl = this->_clusterData->getCluster(_reactant);
	auto dc = cl.getDiffusionCoefficient(gridIndex);
	auto anisotropy =
		this->_clusterData->extraData.anisotropyRatio(_reactant, gridIndex);

	auto clReg = cl.getRegion();
	Composition lo = clReg.getOrigin();

	if (lo.isOnAxis(Species::V)) {
		return dc * 1.0 *
			(::xolotl::core::alphaZrCSinkStrength * anisotropy +
				::xolotl::core::alphaZrASinkStrength /
					(anisotropy * anisotropy));
	}

    //1-D diffusers are assumed to only interact with <a>-type edge dislocation lines
    //The anisotropy factor is assumed equal to 1.0 in this case
	else if (lo.isOnAxis(Species::I)) {
        if (lo[Species::I] < 9) {
            return dc * 1.1 *
                (::xolotl::core::alphaZrCSinkStrength * anisotropy +
                    ::xolotl::core::alphaZrASinkStrength /
                        (anisotropy * anisotropy));
        }
        else if (lo[Species::I] == 9) {
            return dc * 1.1 *
                (::xolotl::core::alphaZrASinkStrength);
        }
    }

	return 1.0;
}
} // namespace network
} // namespace core
} // namespace xolotl
