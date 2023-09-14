#pragma once

#include <xolotl/core/network/impl/SinkReaction.tpp>
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
	const double r1, const double dc0, const double dc1)
{
	constexpr double pi = ::xolotl::core::pi;
	constexpr double rCore = ::xolotl::core::alloyCoreRadius;
	const double zs = 4.0 * pi * (r0 + r1 + rCore);
	using Species = typename TRegion::EnumIndex;
	bool cl0IsSphere = (pairCl0Reg.getOrigin().isOnAxis(Species::V) ||
			 pairCl0Reg.getOrigin().isOnAxis(Species::Void) ||
			 pairCl0Reg.getOrigin().isOnAxis(Species::I)),
		 cl1IsSphere = (pairCl1Reg.getOrigin().isOnAxis(Species::V) ||
			 pairCl1Reg.getOrigin().isOnAxis(Species::Void) ||
			 pairCl1Reg.getOrigin().isOnAxis(Species::I));

	// Simple case
	if (cl0IsSphere && cl1IsSphere) {
		return zs * (dc0 + dc1);
	}

	// Both loops
	double r_l = 0.0, r_s = 0.0;
	if (not cl0IsSphere and not cl1IsSphere) {
		r_l = (r0 < r1) ? r1 : r0;
		r_s = (r0 < r1) ? r0 : r1;
	}
	// Loop and sphere
	else {
		// Which one is sphere?
		r_s = (cl0IsSphere) ? r0 : r1;
		r_l = (cl0IsSphere) ? r1 : r0;
	}

	double p = 1.0 / (1.0 + pow(r_l / (r_s + rCore), 2.0));
	double zl = 4.0 * pi * pi * r_l / log(1.0 + 8.0 * r_l / (r_s + rCore));

	double k_plus = (dc0 + dc1) * (p * zs + (1.0 - p) * zl);
	double bias = 1.0;
	if (pairCl0Reg.getOrigin().isOnAxis(Species::I) ||
		pairCl1Reg.getOrigin().isOnAxis(Species::I)) {
		bias = 1.2;
	}

	return k_plus * bias;
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

	return alloy::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
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

	return alloy::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
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
	if (lo.isOnAxis(Species::Void)) {
		double n = (double)(lo[Species::Void] + hi[Species::Void] - 1) * 0.5;
		//		if (prod1Comp.isOnAxis(Species::I) ||
		// prod2Comp.isOnAxis(Species::I)) { 			be = 3.5 - 3.45 * (pow(n
		// + 1.0, 2.0
		/// 3.0) - pow(n, 2.0 / 3.0));
		//		}
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 1.9 - 3.1 * (cbrt(n * n) - cbrt((n - 1.0) * (n - 1.0)));
		}
	}
	else if (lo.isOnAxis(Species::Faulted)) {
		double n =
			(double)(lo[Species::Faulted] + hi[Species::Faulted] - 1) * 0.5;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 1.9 - 3.2 * (cbrt(n * n) - cbrt((n - 1.0) * (n - 1.0)));
		}
	}
	else if (lo.isOnAxis(Species::V)) {
		double n = (double)(lo[Species::V] + hi[Species::V] - 1) * 0.5;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 1.9 - 3.1 * (cbrt(n * n) - cbrt((n - 1.0) * (n - 1.0)));
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
} // namespace network
} // namespace core
} // namespace xolotl
