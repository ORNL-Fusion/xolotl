#pragma once

#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/core/network/impl/TransformReaction.tpp>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace alloy
{
template <typename TRegion>
KOKKOS_INLINE_FUNCTION
double
getRate(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg, const double r0,
	const double r1, const double dc0, const double dc1, double rdCl[2][2])
{
	constexpr double pi = ::xolotl::core::pi;
	constexpr double rCore = ::xolotl::core::alloyCoreRadius;
	const double zs = 4.0 * pi * (r0 + r1 + rCore);
	const double p = 1.0;

	using Species = typename TRegion::EnumIndex;
	xolotl::core::network::detail::Composition<typename TRegion::VectorType,
		Species>
		lo0 = pairCl0Reg.getOrigin();
	xolotl::core::network::detail::Composition<typename TRegion::VectorType,
		Species>
		lo1 = pairCl1Reg.getOrigin();

	// Determine if clusters are vacancy or interstitial and initialize
	// variables
	bool cl0IsV = lo0.isOnAxis(Species::V) || lo0.isOnAxis(Species::PerfectV) ||
		lo0.isOnAxis(Species::FaultedV);
	bool cl1IsV = lo1.isOnAxis(Species::V) || lo1.isOnAxis(Species::PerfectV) ||
		lo1.isOnAxis(Species::FaultedV);
	double n0 = 0; // size of cluster 0
	double n1 = 0; // size of cluster 1
	double Pl = 1.0; // Capture efficiency for diffusing defect
	double Pli = 1.0; // Capture efficiency for interstitial a-loop
	double Plv = 1.0; // Capture efficiency for vacancy a-loops

	// Determine parameters for cluster 0 based on cluster type and size
	n0 = lo0[Species::V] + lo0[Species::PerfectV] + lo0[Species::FaultedV] +
		lo0[Species::I] + lo0[Species::PerfectI] + lo0[Species::FaultedI];
	n1 = lo1[Species::V] + lo1[Species::PerfectV] + lo1[Species::FaultedV] +
		lo1[Species::I] + lo1[Species::PerfectI] + lo1[Species::FaultedI];

	bool cl0IsSphere = (pairCl0Reg.getOrigin().isOnAxis(Species::V) ||
			 pairCl0Reg.getOrigin().isOnAxis(Species::I)),
		 cl1IsSphere = (pairCl1Reg.getOrigin().isOnAxis(Species::V) ||
			 pairCl1Reg.getOrigin().isOnAxis(Species::I));

	// Cluster 0 is a dislocation loop
	if (not cl0IsSphere) {
		// Define the dislocation capture radius, transition coefficient, and
		// then calculate the reaction rate
		double rd = rdCl[0][cl1IsV];
		double alpha = pow(1 + pow(r0 / (3 * (r1 + rd)), 2), -1);
		double rateSpherical = 4.0 * pi * (r0 + r1 + rd);
		double rateToroidal =
			(4.0 * pi * pi * r0) / log(1 + (8 * r0) / (r1 + rd));

		// Calculate the capture efficiency (assuming only prismatic loops)
		if (cl0IsV)
			Pl = 0.78 * pow(p, -2) + 0.66 * p - 0.44;
		else
			Pl = 0.70 * pow(p, -2) + 0.78 * p - 0.47;

		return ((1 - alpha) * rateToroidal * Pl + alpha * rateSpherical) *
			(dc0 + dc1);
	}

	// Cluster 1 is a dislocation loop:
	if (not cl1IsSphere) {
		// Define the dislocation capture radius, transition coefficient, and
		// then calculate the reaction rate
		double rd = rdCl[1][cl0IsV];
		double alpha = pow(1 + pow(r1 / (3 * (r0 + rd)), 2), -1);
		double rateSpherical = 4.0 * pi * (r0 + r1 + rd);
		double rateToroidal =
			(4.0 * pi * pi * r1) / log(1 + (8 * r1) / (r0 + rd));

		// Calculate the capture efficiency (assuming only prismatic loops)
		if (cl1IsV)
			Pl = 0.78 * pow(p, -2) + 0.66 * p - 0.44;
		else
			Pl = 0.70 * pow(p, -2) + 0.78 * p - 0.47;

		return ((1 - alpha) * rateToroidal * Pl + alpha * rateSpherical) *
			(dc0 + dc1);
	}

	// None of the clusters are loops (interaction is based on spherical volume)
	return zs * (dc0 + dc1);
}
} // namespace alloy

KOKKOS_INLINE_FUNCTION
double
AlloyProductionReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	// Create an array with all possible dislocation capture radii
	// rdCl = {(rdI for cl0, rdV for cl0), (rdI for cl1, rdV for cl1)}
	double rdCl[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
	rdCl[0][0] = this->_clusterData->extraData.dislocationCaptureRadius(
		_reactants[0], 0);
	rdCl[0][1] = this->_clusterData->extraData.dislocationCaptureRadius(
		_reactants[0], 1);
	rdCl[1][0] = this->_clusterData->extraData.dislocationCaptureRadius(
		_reactants[1], 0);
	rdCl[1][1] = this->_clusterData->extraData.dislocationCaptureRadius(
		_reactants[1], 1);

	auto rate = alloy::getRate(
		cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1, rdCl);

	// Divide the rate by 2 for I + I -> faulted or perfect
	if (cl0.getRegion().getOrigin().isOnAxis(Species::I) and
		cl1.getRegion().getOrigin().isOnAxis(Species::I)) {
		auto prod = this->_clusterData->getCluster(_products[0]);
		if (not prod.getRegion().getOrigin().isOnAxis(Species::I)) {
			return rate * 0.5;
		}
	}

	return rate;
}

KOKKOS_INLINE_FUNCTION
double
AlloyDissociationReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	// Create an array with all possible dislocation capture radii
	// rdCl = {(rdI for cl0, rdV for cl0), (rdI for cl1, rdV for cl1)}
	double rdCl[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
	rdCl[0][0] =
		this->_clusterData->extraData.dislocationCaptureRadius(_products[0], 0);
	rdCl[0][1] =
		this->_clusterData->extraData.dislocationCaptureRadius(_products[0], 1);
	rdCl[1][0] =
		this->_clusterData->extraData.dislocationCaptureRadius(_products[1], 0);
	rdCl[1][1] =
		this->_clusterData->extraData.dislocationCaptureRadius(_products[1], 1);

	auto rate = alloy::getRate(
		cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1, rdCl);

	// Divide the rate by 2 for I + I -> faulted or perfect
	if (cl0.getRegion().getOrigin().isOnAxis(Species::I) and
		cl1.getRegion().getOrigin().isOnAxis(Species::I)) {
		auto react = this->_clusterData->getCluster(_reactant);
		if (not react.getRegion().getOrigin().isOnAxis(Species::I)) {
			return rate * 0.5;
		}
	}

	return rate;
}

KOKKOS_INLINE_FUNCTION
double
AlloyDissociationReaction::computeBindingEnergy(double time)
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
	if (lo.isOnAxis(Species::PerfectV)) {
		double n =
			(double)(lo[Species::PerfectV] + hi[Species::PerfectV] - 1) * 0.5;
		//		if (prod1Comp.isOnAxis(Species::I) ||
		// prod2Comp.isOnAxis(Species::I)) { 			be = 3.5 - 3.45 * (pow(n
		// + 1.0, 2.0
		/// 3.0) - pow(n, 2.0 / 3.0));
		//		}
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 1.0 - 3.1 * (cbrt(n * n) - cbrt((n - 1.0) * (n - 1.0)));
		}
	}
	else if (lo.isOnAxis(Species::FaultedV)) {
		double n =
			(double)(lo[Species::FaultedV] + hi[Species::FaultedV] - 1) * 0.5;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 1.0 - 3.2 * (cbrt(n * n) - cbrt((n - 1.0) * (n - 1.0)));
		}
	}
	else if (lo.isOnAxis(Species::V)) {
		double n = (double)(lo[Species::V] + hi[Species::V] - 1) * 0.5;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 1.0 - 3.1 * (cbrt(n * n) - cbrt((n - 1.0) * (n - 1.0)));
		}
	}
	//	else if (lo.isOnAxis(Species::I)) {
	//		double n = (double)(lo[Species::I] + hi[Species::I] - 1) / 2.0;
	//		if (prod1Comp.isOnAxis(Species::I) ||
	// prod2Comp.isOnAxis(Species::I)) { 			be = 3.5 - 2.5 * (pow(n, 2.0
	/// 3.0) - pow(n - 1.0, 2.0 / 3.0));
	//		}
	//	}

	return util::min(5.0, util::max(be, 0.1));
}

KOKKOS_INLINE_FUNCTION
double
AlloySinkReaction::getSinkBias()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	double bias = 1.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);

	auto clReg = cl.getRegion();
	if (clReg.isSimplex()) {
		Composition comp = clReg.getOrigin();
		if (comp.isOnAxis(Species::I)) {
			bias = 1.2;
		}
	}

	return bias;
}

KOKKOS_INLINE_FUNCTION
double
AlloySinkReaction::getSinkStrength()
{
	return ::xolotl::core::alloySinkStrength;
}

KOKKOS_INLINE_FUNCTION
double
AlloyTransformReaction::getSize()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	auto cl = this->_clusterData->getCluster(this->_reactant);

	auto clReg = cl.getRegion();
	Composition lo = clReg.getOrigin();
	Composition hi = clReg.getUpperLimitPoint();

	return (lo[Species::FaultedI] + lo[Species::FaultedV] +
			   hi[Species::FaultedI] + hi[Species::FaultedV] - 2.0) *
		0.5;
}

KOKKOS_INLINE_FUNCTION
double
AlloyTransformReaction::getExponent()
{
	return 1.0;
}

KOKKOS_INLINE_FUNCTION
double
AlloyTransformReaction::getBarrier()
{
	return this->_clusterData->barrierEnergy();
}
} // namespace network
} // namespace core
} // namespace xolotl
