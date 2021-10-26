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
    const double r1, const double dc0, const double dc1, double rdCl[2][2], double p)
{
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

    //Determine if clusters are vacancy or interstitial and initialize variables
    bool cl0IsV = lo0.isOnAxis(Species::V);
    bool cl1IsV = lo1.isOnAxis(Species::V);
    double n0 = 0; //size of cluster 0
    double n1 = 0; //size of cluster 1
    double Pl = 1.0; //Capture efficiency for diffusing defect

    //Determine parameters based on cluster type and size
    if (cl0IsV) n0 = lo0[Species::V];
    else n0 = lo0[Species::I];
    bool cl0IsLoop = (n0 > 9);

    if (cl1IsV) n1 = lo1[Species::V];
    else n1 = lo1[Species::I];
    bool cl1IsLoop = (n1 > 9);

    // These rates are only for 3-D mobile diffusers, but for now assume that it works for 1-D diffusers as well
    // Cluster 0 is a dislocation loop:
    if (cl0IsLoop) {
        //Define the dislocation capture radius, transition coefficient, and then calculate the reaction rate
        double rd = rdCl[0][cl1IsV];
        double alpha = pow(1+pow(r0 / (3*(r1 + rd)), 2), -1);
        double rateSpherical = 4.0 * pi * (r0 + r1 + rd);
        double rateToroidal = (4.0 * pi * pi * r0) / log(1 + (8 * r0) / (r1 + rd));

        //Calculate the capture efficiency (assuming only prismatic loops)
        if (cl0IsV) Pl = 0.78 * pow(p,-2) + 0.66 * p - 0.44;
        else Pl = 0.70 * pow(p,-2) + 0.78 * p - 0.47;

        return ((1 - alpha) * rateToroidal * Pl + alpha * rateSpherical) * (dc0 + dc1);

    }

    // Cluster 1 is a dislocation loop:
    if (cl1IsLoop) {
        //Define the dislocation capture radius, transition coefficient, and then calculate the reaction rate
        double rd = rdCl[1][cl0IsV];
        double alpha = pow(1+pow(r1 / (3*(r0 + rd)), 2), -1);
        double rateSpherical = 4.0 * pi * (r0 + r1 + rd);
        double rateToroidal = (4.0 * pi * pi * r1) / log(1 + (8 * r1) / (r0 + rd));

        //Calculate the capture efficiency (assuming only prismatic loops)
        if (cl1IsV) Pl = 0.78 * pow(p,-2) + 0.66 * p - 0.44;
        else Pl = 0.70 * pow(p,-2) + 0.78 * p - 0.47;

        return ((1 - alpha) * rateToroidal * Pl + alpha * rateSpherical) * (dc0 + dc1);
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

    //Determine which cluster is mobile and retrieve its anisotropy ratio
    double p = 0;
    if (dc0>0) p = this->_clusterData->extraData.anisotropyRatio(_reactants[0], gridIndex);
    else if (dc1>0) p = this->_clusterData->extraData.anisotropyRatio(_reactants[1], gridIndex);

    //Create an array with all possible dislocation capture radii
    //rdCl = {(rdI for cl0, rdV for cl0), (rdI for cl1, rdV for cl1)}
    double rdCl[2][2] = {{0.0,0.0}, {0.0,0.0}};
    rdCl[0][0] = this->_clusterData->extraData.dislocationCaptureRadius(_reactants[0], 0);
    rdCl[0][1] = this->_clusterData->extraData.dislocationCaptureRadius(_reactants[0], 1);
    rdCl[1][0] = this->_clusterData->extraData.dislocationCaptureRadius(_reactants[1], 0);
    rdCl[1][1] = this->_clusterData->extraData.dislocationCaptureRadius(_reactants[1], 1);

    return getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1, rdCl, p);
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

    //Determine which cluster is mobile and retrieve its anisotropy ratio
    double p = 0;
    if (dc0>0) p = this->_clusterData->extraData.anisotropyRatio(_products[0], gridIndex);
    else if (dc1>0) p = this->_clusterData->extraData.anisotropyRatio(_products[1], gridIndex);

    //Create an array with all possible dislocation capture radii
    //rdCl = {(rdI for cl0, rdV for cl0), (rdI for cl1, rdV for cl1)}
    double rdCl[2][2] = {{0.0,0.0}, {0.0,0.0}};
    rdCl[0][0] = this->_clusterData->extraData.dislocationCaptureRadius(_products[0], 0);
    rdCl[0][1] = this->_clusterData->extraData.dislocationCaptureRadius(_products[0], 1);
    rdCl[1][0] = this->_clusterData->extraData.dislocationCaptureRadius(_products[1], 0);
    rdCl[1][1] = this->_clusterData->extraData.dislocationCaptureRadius(_products[1], 1);

    return getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1, rdCl, p);
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
	double dc = cl.getDiffusionCoefficient(gridIndex);
	double anisotropy =
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
