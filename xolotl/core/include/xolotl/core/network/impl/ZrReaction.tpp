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
	// TODO: fix the reaction rate
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

	// Cluster 0 is 1D diffuser
	if (cl0Is1D) {
		return 0.0;
	}
	// Cluster 1 is 1D diffuser
	else if (cl1Is1D) {
		return 0.0;
	}

	// None of them is a 1D diffuser
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

	// TODO: Fix the formulas for V and I

	if (lo.isOnAxis(Species::V)) {
		double n = (double)(lo[Species::V] + hi[Species::V] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 0.0 - 0.0 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
		}
	}
	else if (lo.isOnAxis(Species::I)) {
		double n = (double)(lo[Species::I] + hi[Species::I] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
			be = 0.0 - 0.0 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
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

	// TODO: set the right values in the arrays, verify the formulas

	// Anisotropy ratio
	constexpr Kokkos::Array<double, 6> iAnisotropy = {
		0.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	constexpr Kokkos::Array<double, 10> vAnisotropy = {
		0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

	auto cl = this->_clusterData->getCluster(_reactant);
	double dc = cl.getDiffusionCoefficient(gridIndex);

	auto clReg = cl.getRegion();
	Composition lo = clReg.getOrigin();

	if (lo.isOnAxis(Species::V)) {
		return dc * 1.0 *
			(::xolotl::core::alphaZrASinkStrength *
					vAnisotropy[lo[Species::V]] +
				::xolotl::core::alphaZrCSinkStrength /
					(vAnisotropy[lo[Species::V]] *
						vAnisotropy[lo[Species::V]]));
	}
	else if (lo.isOnAxis(Species::I)) {
		if (lo[Species::I] < iAnisotropy.size())
			return dc * 1.1 *
				(::xolotl::core::alphaZrASinkStrength *
						iAnisotropy[lo[Species::I]] +
					::xolotl::core::alphaZrCSinkStrength /
						(iAnisotropy[lo[Species::I]] *
							iAnisotropy[lo[Species::I]]));
		else
			return dc * 1.1 *
				(::xolotl::core::alphaZrASinkStrength +
					::xolotl::core::alphaZrCSinkStrength);
	}

	return 1.0;
}
} // namespace network
} // namespace core
} // namespace xolotl
