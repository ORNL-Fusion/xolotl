#pragma once

#include <xolotl/core/network/impl/NucleationReaction.tpp>
#include <xolotl/core/network/impl/ReSolutionReaction.tpp>
#include <xolotl/core/network/impl/Reaction.tpp>
#include <xolotl/core/network/impl/SinkReaction.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
KOKKOS_INLINE_FUNCTION
double
NEDissociationReaction::computeBindingEnergy()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	constexpr double xeTrapTable[9] = {
		0.0, 4.31, 2.90, 2.02, 1.09, 0.58, 0.13, -0.25, -0.59};

	double be = 5.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);
	auto prod1 = this->_clusterData->getCluster(this->_products[0]);
	auto prod2 = this->_clusterData->getCluster(this->_products[1]);

	auto clReg = cl.getRegion();
	auto prod1Reg = prod1.getRegion();
	auto prod2Reg = prod2.getRegion();
	if (clReg.isSimplex() && prod1Reg.isSimplex() && prod2Reg.isSimplex()) {
		Composition comp = clReg.getOrigin();
		Composition prod1Comp = prod1Reg.getOrigin();
		Composition prod2Comp = prod2Reg.getOrigin();
		if (comp.isOnAxis(Species::Xe)) {
			if (prod1Comp.isOnAxis(Species::Xe) ||
				prod2Comp.isOnAxis(Species::Xe)) {
				if (comp[Species::Xe] == 2)
					be = 0.5;
				else
					be = 1.0;
			}
			if (prod1Comp.isOnAxis(Species::I) ||
				prod2Comp.isOnAxis(Species::I)) {
				be = xeTrapTable[comp[Species::Xe]];
			}
		}
		else if (comp.isOnAxis(Species::V)) {
			auto size = comp[Species::V];
			be = 1.73 -
				2.59 *
					(pow((double)size, 2.0 / 3.0) -
						pow((double)size - 1.0, 2.0 / 3.0));
		}
		else if (comp.isOnAxis(Species::I)) {
			// Nothing
		}
		else {
			// XeV
			if (prod1Comp.isOnAxis(Species::V) ||
				prod2Comp.isOnAxis(Species::V)) {
				auto amtXe = comp[Species::Xe], amtV = comp[Species::V];
				be = 1.73 -
					2.59 *
						(pow((double)amtV, 2.0 / 3.0) -
							pow((double)amtV - 1.0, 2.0 / 3.0)) +
					2.5 * log(1.0 + ((double)amtXe / (double)amtV));
			}
			if (prod1Comp.isOnAxis(Species::I) ||
				prod2Comp.isOnAxis(Species::I)) {
				auto amtXe = comp[Species::Xe], amtV = comp[Species::V];
				be = 4.88 +
					2.59 *
						(pow((double)amtV, 2.0 / 3.0) -
							pow((double)amtV - 1.0, 2.0 / 3.0)) -
					2.5 * log(1.0 + ((double)amtXe / (double)amtV));
			}
		}
	}
	else {
		Composition lo = clReg.getOrigin();
		Composition hi = clReg.getUpperLimitPoint();
		Composition prod1Comp = prod1Reg.getOrigin();
		Composition prod2Comp = prod2Reg.getOrigin();
		// XeV
		auto amtXe = (double)(lo[Species::Xe] + hi[Species::Xe] - 1) / 2.0;
		auto amtV = (double)(lo[Species::V] + hi[Species::V] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 1.73 -
				2.59 * (pow(amtV, 2.0 / 3.0) - pow(amtV - 1.0, 2.0 / 3.0)) +
				2.5 * log(1.0 + (amtXe / amtV));
		}
		if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
			be = 4.88 +
				2.59 * (pow(amtV, 2.0 / 3.0) - pow(amtV - 1.0, 2.0 / 3.0)) -
				2.5 * log(1.0 + (amtXe / amtV));
		}
	}

	return util::min(5.0, util::max(be, -5.0));
}

KOKKOS_INLINE_FUNCTION
double
NESinkReaction::getSinkBias()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	double bias = 1.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);

	auto clReg = cl.getRegion();
	if (clReg.isSimplex()) {
		Composition comp = clReg.getOrigin();
		if (comp.isOnAxis(Species::I)) {
			bias = 1.05;
		}
	}

	return bias;
}

KOKKOS_INLINE_FUNCTION
double
NESinkReaction::getSinkStrength()
{
	auto cl = this->_clusterData->getCluster(this->_reactant);
	double r = cl.getReactionRadius();
	double latticeParameter = this->_clusterData->latticeParameter();
	double r0 = latticeParameter * 0.5 * sqrt(2.0);
	double rho = 0.0003;
	constexpr double pi = ::xolotl::core::pi;

	double strength = -4.0 * pi * rho / log(pi * rho * (r + r0) * (r + r0));

	return strength;
}
} // namespace network
} // namespace core
} // namespace xolotl
