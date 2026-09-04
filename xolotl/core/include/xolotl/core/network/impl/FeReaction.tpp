#pragma once

#include <xolotl/core/network/FeClusterGenerator.h>
#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/core/network/impl/TrapReaction.tpp>
#include <xolotl/util/MathUtils.h>
#include <iostream> // debugging
#include <fstream>  // debugging



namespace xolotl
{
namespace core
{
namespace network
{
namespace fe
{
template <typename TRegion>
KOKKOS_INLINE_FUNCTION
double
getRate(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg, const double r0,
	const double r1, const double dc0, const double dc1)
{
	constexpr double pi = ::xolotl::core::pi;

	double kPlus =
		4.0 * pi * (r0 + r1 + ::xolotl::core::feCrCoreRadius) * (dc0 + dc1);

	return kPlus;
}
} // namespace fe

KOKKOS_INLINE_FUNCTION
double
FeProductionReaction::getRateForProduction(IndexType gridIndex)
{
	// Standard case
	if (not isLargeBubbleReaction) {
		auto cl0 = this->_clusterData->getCluster(_reactants[0]);
		auto cl1 = this->_clusterData->getCluster(_reactants[1]);

		double r0 = cl0.getReactionRadius();
		double r1 = cl1.getReactionRadius();

		double dc0 = cl0.getDiffusionCoefficient(gridIndex);
		double dc1 = cl1.getDiffusionCoefficient(gridIndex);

		return fe::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
	}



	// SSBM reaction
	const auto dummyRegion = Region(Composition{});

	double r0 = 0.0, r1 = 0.0, dc0 = 0.0, dc1 = 0.0;
	Region cl0Reg = dummyRegion, cl1Reg = dummyRegion;
	auto numClusters = this->_clusterData->numClusters;

	if (this->_reactants[0] >= numClusters) {
		r0 = this->_clusterData->bubbleAvRad();
	}
	else {
		auto cl0 = this->_clusterData->getCluster(this->_reactants[0]);
		r0 = cl0.getReactionRadius();
		dc0 = cl0.getDiffusionCoefficient(gridIndex);
		cl0Reg = cl0.getRegion();
	}

	if (this->_reactants[1] >= numClusters) {
		r1 = this->_clusterData->bubbleAvRad();
	}
	else {
		auto cl1 = this->_clusterData->getCluster(this->_reactants[1]);
		r1 = cl1.getReactionRadius();
		dc1 = cl1.getDiffusionCoefficient(gridIndex);
		cl1Reg = cl1.getRegion();
	}
	
	return fe::getRate(cl0Reg, cl1Reg, r0, r1, dc0, dc1);

}

KOKKOS_INLINE_FUNCTION
double
FeDissociationReaction::getRateForProduction(IndexType gridIndex)
{
	// Standard case
	if (not isLargeBubbleReaction) {
		auto cl0 = this->_clusterData->getCluster(this->_products[0]);
		auto cl1 = this->_clusterData->getCluster(this->_products[1]);

		double r0 = cl0.getReactionRadius();
		double r1 = cl1.getReactionRadius();

		double dc0 = cl0.getDiffusionCoefficient(gridIndex);
		double dc1 = cl1.getDiffusionCoefficient(gridIndex);

		return fe::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
	}

	// SSBM reaction
	const auto dummyRegion = Region(Composition{});

	double r0 = 0.0, r1 = 0.0, dc0 = 0.0, dc1 = 0.0;
	Region cl0Reg = dummyRegion, cl1Reg = dummyRegion;
	auto numClusters = this->_clusterData->numClusters;

	if (this->_products[0] >= numClusters) {
		r0 = this->_clusterData->bubbleAvRad();
	}
	else {
		auto cl0 = this->_clusterData->getCluster(this->_products[0]);
		r0 = cl0.getReactionRadius();
		dc0 = cl0.getDiffusionCoefficient(gridIndex);
		cl0Reg = cl0.getRegion();
	}

	if (this->_products[1] >= numClusters) {
		r1 = this->_clusterData->bubbleAvRad();
	}
	else {
		auto cl1 = this->_clusterData->getCluster(this->_products[1]);
		r1 = cl1.getReactionRadius();
		dc1 = cl1.getDiffusionCoefficient(gridIndex);
		cl1Reg = cl1.getRegion();
	}

	return fe::getRate(cl0Reg, cl1Reg, r0, r1, dc0, dc1);

}


KOKKOS_INLINE_FUNCTION
void
FeProductionReaction::computeCoefficients()
{
	// Check if the large bubble is involved
	if (isLargeBubbleReaction) {
		constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
		for (auto i : speciesRangeNoI) {
			this->_widths(i()) = 1.0;
		}
		this->_coefs(0, 0, 0, 0) = 1.0;
	}
	else {
		// Standard case
		return Superclass::computeCoefficients();
	}
}

KOKKOS_INLINE_FUNCTION
void
FeProductionReaction::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	// Standard case
	if (not isLargeBubbleReaction) {
		return Superclass::computeFlux(concentrations, fluxes, gridIndex);
	}

	// The rate need to be computed each time because it depends on the current
	// large bubble size
	auto rate = getRateForProduction(gridIndex);

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
	auto numClusters = this->_clusterData->numClusters;

	// Large bubble is one of the reactants
	if (this->_reactants[0] >= numClusters or
		this->_reactants[1] >= numClusters) {
		// Get the standard cluster
		auto stdClusterId = (this->_reactants[0] >= numClusters) ?
			this->_reactants[1] :
			this->_reactants[0];
		
		auto cl = this->_clusterData->getCluster(stdClusterId);
		auto clReg = cl.getRegion();
		auto orig = clReg.getOrigin();
		Composition comp(orig);
		// Get the SSBM cluster
		auto ssbmId = (this->_reactants[0] >= numClusters) ?
			this->_reactants[0] :
			this->_reactants[1];

		// Compute the flux
		double f = this->_coefs(0, 0, 0, 0) * concentrations(stdClusterId) *
			concentrations(ssbmId) * rate;

		// Vacancy case
		// V_k + B -> B
		if (comp[Species::V] > 0) {
			// The standard cluster always loses the flux
			Kokkos::atomic_sub(&fluxes[stdClusterId], f);

			// The V size increases (kind of cheating on the ID)
			Kokkos::atomic_add(&fluxes[ssbmId + 2], f * comp[Species::V]);
		}

		// Interstitial case
		if (comp[Species::I] > 0) {
			// I_k + B -> B
			// The reaction can happen only if the average V is larger than the
			// threshold and if the H/V ratio is realistic
			// Special factor to deal with the threshold size
				double gauss = 1.0;

				// Special case where the product is not the single size
				if (this->_products[0] < numClusters) {
					// Only if the large void has a specific size
					auto avVoid =
						concentrations(ssbmId + 2) / concentrations(ssbmId);
					if (concentrations(ssbmId) == 0.0)
						avVoid = 0.0;

					// Get the product composition
					auto pr =
						this->_clusterData->getCluster(this->_products[0]);
					auto prReg = pr.getRegion();
					Composition prComp(prReg.getOrigin());

					// Target value for the reaction to happen
					double target = prComp[Species::V] + comp[Species::I];

					// Gaussian function around it
					double twoSigmaTwo = 0.1;
					gauss = exp(-(avVoid - target) * (avVoid - target) /
								twoSigmaTwo) /
						sqrt(::xolotl::core::pi * twoSigmaTwo);

					// The large bubble concentration decreases
					Kokkos::atomic_sub(&fluxes[ssbmId], f * gauss);
					// The product concentration increases
					Kokkos::atomic_add(&fluxes[this->_products[0]], f * gauss);

					// The V size decreases even more
					Kokkos::atomic_sub(&fluxes[ssbmId + 2],
						f *
							(prComp[Species::V]) *
							gauss);
				}

				// In every case

				// The standard cluster always loses the flux
				Kokkos::atomic_sub(&fluxes[stdClusterId], f * gauss);

				// The V size decreases
				Kokkos::atomic_sub(
					&fluxes[ssbmId + 2], f * comp[Species::I] * gauss);
			}
		

		// H case
		// H_k + B -> B
			if (comp[Species::He] > 0) {
				// The standard cluster always loses the flux
			Kokkos::atomic_sub(&fluxes[stdClusterId], f);

			// The He size increases
			Kokkos::atomic_add(&fluxes[ssbmId + 1], f * comp[Species::He]);
			}
	}

	// Large bubble is one of the product
	else {
		auto cR1 = concentrations[this->_reactants[0]];
		auto cR2 = concentrations[this->_reactants[1]];
		auto cl1 = this->_clusterData->getCluster(this->_reactants[0]);
		auto cl1Reg = cl1.getRegion();
		auto orig1 = cl1Reg.getOrigin();
		Composition comp1(orig1);
		auto cl2 = this->_clusterData->getCluster(this->_reactants[1]);
		auto cl2Reg = cl2.getRegion();
		auto orig2 = cl2Reg.getOrigin();
		Composition comp2(orig2);

		double f = this->_coefs(0, 0, 0, 0) * cR1 * cR2 * rate;

		// Vacancy case
		// V_a + HV_b -> B
		if (orig1.isOnAxis(Species::V) or orig2.isOnAxis(Species::V)) {
			// Compute the total size
			auto totalSize = comp1[Species::V] + comp2[Species::V];
			// Both reactants decrease
			Kokkos::atomic_sub(&fluxes[this->_reactants[0]], f);
			Kokkos::atomic_sub(&fluxes[this->_reactants[1]], f);
			// The large bubble increases, as well as average V and H
			Kokkos::atomic_add(&fluxes[this->_products[0]], f);
			Kokkos::atomic_add(&fluxes[this->_products[0] + 2], f * totalSize);
				totalSize = comp1[Species::He] + comp2[Species::He];
				Kokkos::atomic_add(
					&fluxes[this->_products[0] + 1], f * totalSize);
		}

		// H case
		// H_a + H_bV -> B
			if (orig1.isOnAxis(Species::He) or orig2.isOnAxis(Species::He)) {
				// Compute the total size
				auto totalDSize = comp1[Species::He] + comp2[Species::He];
				auto totalVSize = comp1[Species::V] + comp2[Species::V];
				// Both reactants decrease
				Kokkos::atomic_sub(&fluxes[this->_reactants[0]], f);
				Kokkos::atomic_sub(&fluxes[this->_reactants[1]], f);
				// The large bubble increases, as well as average H and V
				Kokkos::atomic_add(&fluxes[this->_products[0]], f);
				Kokkos::atomic_add(
					&fluxes[this->_products[0] + 1], f * totalDSize);
				Kokkos::atomic_add(
					&fluxes[this->_products[0] + 2], f * totalVSize);
			}
	}
}

KOKKOS_INLINE_FUNCTION
void
FeProductionReaction::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	// Standard case
	if (not isLargeBubbleReaction) {
		return Superclass::computePartialDerivatives(
			concentrations, values, gridIndex);
	}

	// The _connEntries for large bubbles are ordered in the same way as the
	// moments He H V
	// The rate need to be computed each time because it depends on the current
	// large bubble size
	auto rate = getRateForProduction(gridIndex);

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
	auto numClusters = this->_clusterData->numClusters;

	// Large bubble is one of the reactants
	if (this->_reactants[0] >= numClusters or
		this->_reactants[1] >= numClusters) {
		// Get the standard cluster
		auto stdClusterId = (this->_reactants[0] >= numClusters) ?
			this->_reactants[1] :
			this->_reactants[0];
		auto cl = this->_clusterData->getCluster(stdClusterId);
		auto clReg = cl.getRegion();
		auto orig = clReg.getOrigin();
		Composition comp(orig);
		// Get the SSBM cluster
		auto ssbmId = (this->_reactants[0] >= numClusters) ?
			this->_reactants[0] :
			this->_reactants[1];
		// Get the concentrations
		auto stdC = concentrations(stdClusterId);
		auto bC = concentrations(ssbmId);

		// Compute the flux
		double f = this->_coefs(0, 0, 0, 0) * rate;

		// Vacancy case
		if (comp[Species::V] > 0) {
			// V_k + B -> B
			// The standard cluster always loses the flux
			if (this->_reactants[0] >= numClusters) {
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][0][0][0]), f * stdC);
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][0][1][0]), f * bC);
			}
			else {
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][0][1][0]), f * stdC);
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][0][0][0]), f * bC);
			}

			// The V size increases
			f = this->_coefs(0, 0, 0, 0) * rate * comp[Species::V];
			if (this->_reactants[0] >= numClusters) {
				Kokkos::atomic_add(
					&values(this->_connEntries[0][3][0][0]), f * stdC);
				Kokkos::atomic_add(
					&values(this->_connEntries[0][3][1][0]), f * bC);
			}
			else {
				Kokkos::atomic_add(
					&values(this->_connEntries[1][3][1][0]), f * stdC);
				Kokkos::atomic_add(
					&values(this->_connEntries[1][3][0][0]), f * bC);
			}

			//_connEntries[][][][] understanding still shaky, looked over 800H, H-Blister, and Reaction.cpp files
			//the work shown here is the copied directly from H-Blister file, although I think we may need to change this
			//if instability still occur, this is my first suspect
			// current understanding of _connEntries[a][b][c][d]
			// a is the term we are looking at in the reaction ~ A+B -> C
			// A is 0, B is 1, C is 2
			// b is the derivative we want to take respect to
			// Trap is 0, He is 1, V is 2, I is 3 (in H-Blister different species are involved thus a differnt ordering)
			// c is the moment ordering 
			// 0 is zeroth moment, 1 is first moment, 2 is second moment
			// d is something I am not sure of
			// Please correct if I am missing something/everything
		}
		// Interstitial case
		if (comp[Species::I] > 0) {
			// I_k + B -> B
			// Special factor to deal with the threshold size
				double gauss = 1.0;
			
			// The reaction can happen only if the average V is larger than the
			// threshold and if the H/V ratio is realistic
			auto avV = concentrations(ssbmId + 2) / concentrations(ssbmId);
			auto avH = concentrations(ssbmId + 1) / concentrations(ssbmId);
			if (concentrations(ssbmId) == 0.0) {
				avV = 0.0;
				avH = 0.0;
			}

			// Get the product composition
					auto pr =
						this->_clusterData->getCluster(this->_products[0]);
					auto prReg = pr.getRegion();
					Composition prComp(prReg.getOrigin());

					// Target value for the reaction to happen
					double target = prComp[Species::V] + comp[Species::I];

					// Gaussian function around it
					double twoSigmaTwo = 0.1;
					gauss = exp(-(avV - target) * (avV - target) /
									   twoSigmaTwo) /
						sqrt(::xolotl::core::pi * twoSigmaTwo);


			// The standard cluster always loses the flux
			f = this->_coefs(0, 0, 0, 0) * rate * gauss;
			if (this->_reactants[0] >= numClusters) {
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][0][0][0]), f * stdC);
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][0][1][0]), f * bC);
			}
			else {
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][0][1][0]), f * stdC);
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][0][0][0]), f * bC);
			}

			// The V size decreases
			f = this->_coefs(0, 0, 0, 0) * rate * comp[Species::I] * gauss;
			if (this->_reactants[0] >= numClusters) {
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][3][0][0]), f * stdC);
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][3][1][0]), f * bC);
			}
			else {
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][3][1][0]), f * stdC);
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][3][0][0]), f * bC);
			}
		}

		// He case
			if (comp[Species::He] > 0) {
				// He_k + B -> B

				// The standard cluster always loses the flux
				if (this->_reactants[0] >= numClusters) {
					Kokkos::atomic_sub(
						&values(this->_connEntries[1][0][0][0]), f * stdC);
					Kokkos::atomic_sub(
						&values(this->_connEntries[1][0][1][0]), f * bC);
				}
				else {
					Kokkos::atomic_sub(
						&values(this->_connEntries[0][0][1][0]), f * stdC);
					Kokkos::atomic_sub(
						&values(this->_connEntries[0][0][0][0]), f * bC);
				}

				// The He size increases
				f = this->_coefs(0, 0, 0, 0) * rate * comp[Species::He];
				if (this->_reactants[0] >= numClusters) {
					Kokkos::atomic_add(
						&values(this->_connEntries[0][2][0][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[0][2][1][0]), f * bC);
				}
				else {
					Kokkos::atomic_add(
						&values(this->_connEntries[1][2][1][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[1][2][0][0]), f * bC);
				}

			}
	}

	// Large bubble is one of the product
	else {
		auto cR1 = concentrations[this->_reactants[0]];
		auto cR2 = concentrations[this->_reactants[1]];
		auto cl1 = this->_clusterData->getCluster(this->_reactants[0]);
		auto cl1Reg = cl1.getRegion();
		auto orig1 = cl1Reg.getOrigin();
		Composition comp1(orig1);
		auto cl2 = this->_clusterData->getCluster(this->_reactants[1]);
		auto cl2Reg = cl2.getRegion();
		auto orig2 = cl2Reg.getOrigin();
		Composition comp2(orig2);

		double f = this->_coefs(0, 0, 0, 0) * rate;

		// Vacancy case
		// V_a + H_bV -> B
		if (orig1.isOnAxis(Species::V) or orig2.isOnAxis(Species::V)) {
			// Both reactants decrease
			Kokkos::atomic_sub(
				&values(this->_connEntries[0][0][0][0]), f * cR2);
			Kokkos::atomic_sub(
				&values(this->_connEntries[1][0][0][0]), f * cR2);
			Kokkos::atomic_sub(
				&values(this->_connEntries[0][0][1][0]), f * cR1);
			Kokkos::atomic_sub(
				&values(this->_connEntries[1][0][1][0]), f * cR1);

			// The large bubble increases, as well as average H and V
			Kokkos::atomic_add(
				&values(this->_connEntries[2][0][0][0]), f * cR2);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][0][1][0]), f * cR1);
				f = this->_coefs(0, 0, 0, 0) * rate *
					(comp1[Species::He] + comp2[Species::He]);
				Kokkos::atomic_add(
					&values(this->_connEntries[2][2][0][0]), f * cR2);
				Kokkos::atomic_add(
					&values(this->_connEntries[2][2][1][0]), f * cR1);
			f = this->_coefs(0, 0, 0, 0) * rate *
				(comp1[Species::V] + comp2[Species::V]);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][3][0][0]), f * cR2);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][3][1][0]), f * cR1);
		}

		// He case
		// He_a + He_bV -> B
			if (orig1.isOnAxis(Species::He) or orig2.isOnAxis(Species::He)) {
				// Both reactants decrease
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][0][0][0]), f * cR2);
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][0][0][0]), f * cR2);
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][0][1][0]), f * cR1);
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][0][1][0]), f * cR1);

				// The large bubble increases, as well as average H and V
				Kokkos::atomic_add(
					&values(this->_connEntries[2][0][0][0]), f * cR2);
				Kokkos::atomic_add(
					&values(this->_connEntries[2][0][1][0]), f * cR1);
				f = this->_coefs(0, 0, 0, 0) * rate *
					(comp1[Species::He] + comp2[Species::He]);
				Kokkos::atomic_add(
					&values(this->_connEntries[2][2][0][0]), f * cR2);
				Kokkos::atomic_add(
					&values(this->_connEntries[2][2][1][0]), f * cR1);
				f = this->_coefs(0, 0, 0, 0) * rate *
					(comp1[Species::V] + comp2[Species::V]);
				Kokkos::atomic_add(
					&values(this->_connEntries[2][3][0][0]), f * cR2);
				Kokkos::atomic_add(
					&values(this->_connEntries[2][3][1][0]), f * cR1);
			}
	}
}

KOKKOS_INLINE_FUNCTION
void
FeDissociationReaction::computeCoefficients()
{
	// Check if the large bubble is involved
	if (isLargeBubbleReaction) {
		constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
		for (auto i : speciesRangeNoI) {
			this->_widths(i()) = 1.0;
		}
		this->_coefs(0, 0, 0, 0) = 1.0;
	}
	else {
		// Standard case
		Superclass::computeCoefficients();
	}
}

KOKKOS_INLINE_FUNCTION
void
FeDissociationReaction::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	// Standard case
	if (not isLargeBubbleReaction) {
		return Superclass::computeFlux(concentrations, fluxes, gridIndex);
	}

	// The rate need to be computed each time because it depends on the current
	// large bubble size
	auto rate = getRateForProduction(gridIndex);

	// Binding energy
	double omega = this->_clusterData->atomicVolume();
	double T = this->_clusterData->temperature(gridIndex);
	constexpr double k_B = ::xolotl::core::kBoltzmann;
	double E_b = 1.9; // He in Void converged binding energy value

	rate *= (1.0 / omega) * std::exp(-E_b / (k_B * T));

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
	auto numClusters = this->_clusterData->numClusters;

	// Get the standard cluster
	auto stdClusterId = (this->_products[0] >= numClusters) ?
		this->_products[1] :
		this->_products[0];
	auto cl = this->_clusterData->getCluster(stdClusterId);
	auto clReg = cl.getRegion();
	auto orig = clReg.getOrigin();
	Composition comp(orig);

	// Get the SSBM cluster
	auto ssbmId = (this->_products[0] >= numClusters) ? this->_products[0] : this->_products[1];

	// Compute the flux
	double f = this->_coefs(0, 0, 0, 0) * rate * concentrations(stdClusterId) * concentrations(ssbmId);

	// Large bubble is the reactant
	if (this->_reactant >= numClusters) {
		// H case
		// B -> H_1 + B
			if (comp[Species::He] > 0) {
				// The standard cluster always gains the flux
				Kokkos::atomic_add(&fluxes[stdClusterId], f);

				// The average H decrease
				Kokkos::atomic_sub(&fluxes[ssbmId + 1], f * comp[Species::He]);
			}
	}
}

KOKKOS_INLINE_FUNCTION
void
FeDissociationReaction::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	// Standard case
	if (not isLargeBubbleReaction) {
		return Superclass::computePartialDerivatives(
			concentrations, values, gridIndex);
	}

	// The rate need to be computed each time because it depends on the current
	// large bubble size
	auto rate = getRateForProduction(gridIndex);

	// Binding energy
	double omega = this->_clusterData->atomicVolume();
	double T = this->_clusterData->temperature(gridIndex);
	constexpr double k_B = ::xolotl::core::kBoltzmann;
	double E_b = 1.9; // H in Void

	rate *= (1.0 / omega) * std::exp(-E_b / (k_B * T));

	// Compute the flux
	double f = this->_coefs(0, 0, 0, 0) * rate;

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
	auto numClusters = this->_clusterData->numClusters;

	// Get the standard cluster
	auto stdClusterId = (this->_products[0] >= numClusters) ?
		this->_products[1] :
		this->_products[0];
	auto cl = this->_clusterData->getCluster(stdClusterId);
	auto clReg = cl.getRegion();
	auto orig = clReg.getOrigin();
	Composition comp(orig);

	// Get the SSBM cluster
	auto ssbmId = (this->_products[0] >= numClusters) ? this->_products[0] : this->_products[1];

	// Get the concentration
	auto bC = concentrations(ssbmId);

	// Large bubble is the reactant
	if (this->_reactant >= numClusters) {
		// H case
		// B -> H_1 + B
			if (comp[Species::He] > 0) {
				// The standard cluster always gains the flux
				if (this->_products[0] >= numClusters) {
					Kokkos::atomic_add(
						&values(this->_connEntries[2][0][0][0]), f * bC);
				}
				else {
					Kokkos::atomic_add(
						&values(this->_connEntries[1][0][0][0]), f * bC);
				}

				// The H size decreases
				f = this->_coefs(0, 0, 0, 0) * rate * comp[Species::He];
				if (this->_products[0] >= numClusters) {
					Kokkos::atomic_sub(
						&values(this->_connEntries[1][2][0][0]), f * bC);
				}
				else {
					Kokkos::atomic_sub(
						&values(this->_connEntries[2][2][0][0]), f * bC);
				}
			}
	}
}

KOKKOS_INLINE_FUNCTION
double
FeDissociationReaction::computeBindingEnergy(double time)
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	constexpr double heTrapTable[9] = {
		0.0, 4.31, 2.90, 2.02, 1.09, 0.58, 0.13, -0.25, -0.59};

	constexpr double heHeTable[9] = {
		0.0, 0.0, 0.43, 0.95, 0.98, 1.0, 1.0, 1.0, 1.0};

	constexpr double vVTable[5] = {0.0, 0.0, 0.30, 0.37, 0.62};

	constexpr double bubbleHeTable[4][4] = {{2.3, 2.84, 3.3, 3.84},
		{1.84, 2.75, 2.96, 3.12}, {1.83, 2.07, 2.91, 3.16},
		{1.91, 2.36, 2.57, 3.05}};

	constexpr double bubbleV1Table[8] = {
		3.71, 4.59, 5.52, 6.02, 6.48, 6.86, 7.20, 7.49};

	constexpr double bubbleVTable[3][4] = {{0.78, 1.61, 1.85, 2.30},
		{0.83, 1.04, 1.80, 2.03}, {1.16, 1.32, 1.57, 1.97}};

	constexpr double bubbleITable[3][4] = {{5.83, 5.0, 4.76, 4.31},
		{5.78, 5.57, 4.81, 4.58}, {5.45, 5.29, 5.04, 4.64}};

	double be = 5.0;

//	std::ofstream outFile("FeReaction.logOut",std::ios::app);
//        outFile << "In Binding Energy :\n";

	
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
		if (comp.isOnAxis(Species::He)) {
			if (prod1Comp.isOnAxis(Species::He) ||
				prod2Comp.isOnAxis(Species::He)) {
				be = heHeTable[comp[Species::He]];
			}
			if (prod1Comp.isOnAxis(Species::I) ||
				prod2Comp.isOnAxis(Species::I)) {
				be = heTrapTable[comp[Species::He]];
			}
		}
		else if (comp.isOnAxis(Species::V)) {
			auto size = comp[Species::V];
			if (size < 5)
				be = vVTable[size];
			else {
				be = 1.73 -
					2.59 *
						(pow((double)size, 2.0 / 3.0) -
							pow((double)size - 1.0, 2.0 / 3.0));
			}
		}
		else if (comp.isOnAxis(Species::I)) {
			// Nothing
		}
		else {
			// HeV
			auto amtHe = comp[Species::He];
			auto amtV = comp[Species::V];
			// HeV -> V
			if (prod1Comp.isOnAxis(Species::V) ||
				prod2Comp.isOnAxis(Species::V)) {

				if (amtV == 1 and amtHe < 9)
					be = bubbleV1Table[amtHe - 1];
				else if (amtV < 5 and amtHe < 5) {
					be = bubbleVTable[amtV - 2][amtHe - 1];
				}
				else
				{
					be = 1.73 -
						2.59 *
							(pow((double)amtV, 2.0 / 3.0) -
								pow((double)amtV - 1.0, 2.0 / 3.0)) +
						2.50 * log( 1 + ((double)amtHe) /((double)amtV));
				}
			}
			// HeV -> I
			if (prod1Comp.isOnAxis(Species::I) ||
				prod2Comp.isOnAxis(Species::I)) {
				if (amtV < 4 and amtHe < 5)
					be = bubbleITable[amtV - 1][amtHe - 1];
				else
				{
					be = 4.88 +
						2.59 *
							(pow((double)amtV, 2.0 / 3.0) -
								pow((double)amtV - 1.0, 2.0 / 3.0)) -
						2.50 * log( 1 + ((double)amtHe) /((double)amtV));
				}
			}
		}
	}
	else {
		Composition lo = clReg.getOrigin();
		Composition hi = clReg.getUpperLimitPoint();
		Composition prod1Comp = prod1Reg.getOrigin();
		Composition prod2Comp = prod2Reg.getOrigin();
		// HeV
		auto amtHe = (double)(lo[Species::He] + hi[Species::He] - 1) / 2.0;
		auto amtV = (double)(lo[Species::V] + hi[Species::V] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 1.73 -
				2.59 * (pow(amtV, 2.0 / 3.0) - pow(amtV - 1.0, 2.0 / 3.0)) +
				2.50 * log( 1 + ((double)amtHe) /((double)amtV));
		}
		if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
			be = 4.88 +
				2.59 * (pow(amtV, 2.0 / 3.0) - pow(amtV - 1.0, 2.0 / 3.0)) -
				2.50 * log( 1 + ((double)amtHe) /((double)amtV));
		}
	}

	return util::min(5.0, util::max(be, -5.0));
	
}

KOKKOS_INLINE_FUNCTION
double
FeSinkReaction::getSinkBias()
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
		if (comp.isOnAxis(Species::He)) {
			bias = 0.8;
		}
	}

	return bias;
}

KOKKOS_INLINE_FUNCTION
double
FeSinkReaction::getSinkStrength()
{
	auto cl = this->_clusterData->getCluster(this->_reactant);
	double r = cl.getReactionRadius();
	double latticeParameter = this->_clusterData->latticeParameter();
	double r0 = latticeParameter * 0.75 * sqrt(3.0);
	double rho = 0.00025;
	constexpr double pi = ::xolotl::core::pi;

	double strength =
		-4.0 * pi * rho * (r + r0) / log(pi * rho * (r + r0) * (r + r0));

	return strength;
}

KOKKOS_INLINE_FUNCTION
IdType
FeTrapReaction::getId()
{
	using Composition = typename Superclass::Composition;

	// Here the id is the size of the trap
	auto cl = this->_clusterData->getCluster(this->_trapped);
	Composition comp = cl.getRegion().getOrigin();
	return comp[Species::Trap] - 1;
}

KOKKOS_INLINE_FUNCTION
double
FeTrapReaction::getStrength()
{
	return strength;
}

KOKKOS_INLINE_FUNCTION
double
FeTrapReaction::getEnergy()
{
	return energy;
}

KOKKOS_INLINE_FUNCTION
double
FeTrapReaction::getFrequency()
{
	return frequency;
}

KOKKOS_INLINE_FUNCTION
double
FeTrapReaction::getDensity()
{
	return density;
}
} // namespace network
} // namespace core
} // namespace xolotl
