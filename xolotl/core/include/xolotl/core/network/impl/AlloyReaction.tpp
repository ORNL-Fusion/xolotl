#pragma once

#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
KOKKOS_INLINE_FUNCTION
double
AlloyDissociationReaction::computeBindingEnergy()
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
		double n = (double)(lo[Species::Void] + hi[Species::Void] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
			be = 3.5 - 3.45 * (pow(n + 1.0, 2.0 / 3.0) - pow(n, 2.0 / 3.0));
		}
		else if (prod1Comp.isOnAxis(Species::V) ||
			prod2Comp.isOnAxis(Species::V)) {
			be = 1.9 - 3.1 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
		}
	}
	else if (lo.isOnAxis(Species::Faulted)) {
		double n =
			(double)(lo[Species::Faulted] + hi[Species::Faulted] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 1.9 - 3.2 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
		}
	}
	else if (lo.isOnAxis(Species::V)) {
		double n = (double)(lo[Species::V] + hi[Species::V] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 1.9 - 3.1 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
		}
	}
	else if (lo.isOnAxis(Species::I)) {
		double n = (double)(lo[Species::I] + hi[Species::I] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
			be = 3.5 - 2.5 * (pow(n, 2.0 / 3.0) - pow(n - 1.0, 2.0 / 3.0));
		}
	}

	return util::max(0.1, be);
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
	return ::xolotl::core::alloysinkStrength;
}
} // namespace network
} // namespace core
} // namespace xolotl
