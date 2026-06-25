#pragma once

#include <xolotl/core/network/FeClusterGenerator.h>
#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/core/network/impl/TrapReaction.tpp>
#include <xolotl/util/MathUtils.h>

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
	auto cl0 = this->_clusterData->getCluster(_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return fe::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
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
		Superclass::computeCoefficients();
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
		// Get the index to switch between V and I type SSBM
		auto typeId = (ssbmId - numClusters) / 2;

		// Compute the flux
		double f = this->_coefs(0, 0, 0, 0) * concentrations(stdClusterId) *
			concentrations(ssbmId) * rate;

		// Vacancy case
		// V_k + B -> B
		if (comp[Species::V] > 0) {
			// SSBM V case
				// The standard cluster always loses the flux
				Kokkos::atomic_sub(&fluxes[stdClusterId], f);

				// The V size increases (kind of cheating on the ID)
				Kokkos::atomic_add(&fluxes[ssbmId + 1], f * comp[Species::V]);
		}

		// Interstitial case
		if (comp[Species::I] > 0) {
				// Special factor to deal with the threshold size
				double gauss = 1.0;

				// Special case where the product is not the single size
				if (this->_products[0] < numClusters) {
					// I_a + B -> HV_b
					// Only if the large void has a specific size
					auto avVoid =
						concentrations(ssbmId + 1) / concentrations(ssbmId);
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
					Kokkos::atomic_sub(&fluxes[ssbmId + 1],
						f * prComp[Species::V] * gauss);
				}

				// In every case

				// The standard cluster always loses the flux
				Kokkos::atomic_sub(&fluxes[stdClusterId], f * gauss);

				// The V size decreases
				Kokkos::atomic_sub(
					&fluxes[ssbmId + 1], f * comp[Species::I] * gauss);
		}

		// Helium case
		if (comp[Species::He] > 0) {
			// The standard cluster always loses the flux
			Kokkos::atomic_sub(&fluxes[stdClusterId], f);

			// The He size increases
			Kokkos::atomic_add(&fluxes[ssbmId + 2], f * comp[Species::He]);
/*
				// Trap mutation case
				// Compute the average concentrations
				auto conc = concentrations(ssbmId);
				auto avV = concentrations(ssbmId + 1) / conc;
				auto avHe = concentrations(ssbmId + 2) / conc;
				if (conc == 0.0) {
					avV = 0.0;
					avHe = 0.0;
				}

				// Sigmoid
				double x = ((comp[Species::He] + avHe) /
							   fe::getMaxHePerV(avV,avHe/avV)) -
					1.0;
				double sigmo = 1.0 / (1.0 + std::exp(-50.0 * x));
				double nI = 1.0;

				// The other product increases
				if (this->_products[1] != Superclass::invalidIndex) {
					Kokkos::atomic_add(
						&fluxes[this->_products[1]], nI * f * sigmo);

					// The average V increases
					Composition prodComp(
						this->_clusterData->getCluster(this->_products[1])
							.getRegion()
							.getOrigin());
					Kokkos::atomic_add(&fluxes[ssbmId + 1],
						nI * f * prodComp[Species::I] * sigmo);
				}
				else {
					Kokkos::atomic_add(&fluxes[ssbmId + 1], nI * f * sigmo);
				}
*/
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
			// The large bubble increases, as well as average V
			Kokkos::atomic_add(&fluxes[this->_products[0]], f);
			Kokkos::atomic_add(&fluxes[this->_products[0] + 1], f * totalSize);
		}

		// Interstitial case
		if (orig1.isOnAxis(Species::I) or orig2.isOnAxis(Species::I)) {
			// Compute the total size
			auto totalSize = comp1[Species::I] + comp2[Species::I];
			// Both reactants decrease
			Kokkos::atomic_sub(&fluxes[this->_reactants[0]], f);
			Kokkos::atomic_sub(&fluxes[this->_reactants[1]], f);
			// The large bubble increases, as well as average I
			Kokkos::atomic_add(&fluxes[this->_products[0]], f);
			Kokkos::atomic_add(&fluxes[this->_products[0] + 1], f * totalSize);
		}

		// Helium case
		// He_a + He_bV -> B
		if (orig1.isOnAxis(Species::He) or orig2.isOnAxis(Species::He)) {
			// Compute the total size
			auto totalHeSize = comp1[Species::He] + comp2[Species::He];
			auto totalVSize = comp1[Species::V] + comp2[Species::V];
			// Both reactants decrease
			Kokkos::atomic_sub(&fluxes[this->_reactants[0]], f);
			Kokkos::atomic_sub(&fluxes[this->_reactants[1]], f);
			// The large bubble increases, as well as average He and V
			Kokkos::atomic_add(&fluxes[this->_products[0]], f);
			Kokkos::atomic_add(&fluxes[this->_products[0] + 2], f * totalHeSize);
			Kokkos::atomic_add(&fluxes[this->_products[0] + 1], f * totalVSize);
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
		// Get the index to switch between V and I type SSBM
		auto typeId = (ssbmId - numClusters) / 2;

		// Compute the flux
		double f = this->_coefs(0, 0, 0, 0) * rate;

		// Vacancy case
		if (comp[Species::V] > 0) {
			// SSBM V case
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
						&values(this->_connEntries[0][1][0][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[0][1][1][0]), f * bC);
				}
				else {
					Kokkos::atomic_add(
						&values(this->_connEntries[1][1][1][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[1][1][0][0]), f * bC);
				}
		}
		// Interstitial case
/*		if (comp[Species::I] > 0) {
			// SSBM V case
				// Special factor to deal with the threshold size
				double gauss = 1.0;

				// Special case where the product is not the single size
				if (this->_products[0] < numClusters) {
					// Only if the large void has a specific size
					auto avVoid =
						concentrations(ssbmId + 1) / concentrations(ssbmId);
					if (concentrations(ssbmId) == 0.0)
						avVoid = 0.0;

					// Get the product composition
					auto pr =
						this->_clusterData->getCluster(this->_products[0]);
					auto prReg = pr.getRegion();
					Composition prComp(prReg.getOrigin());

					// Target value for the reaction to happen
					double target = prComp[Species::V] +
						comp[Species::I];

					// Gaussian function around it
					double twoSigmaTwo = 0.1;
					double gauss = exp(-(avVoid - target) * (avVoid - target) /
									   twoSigmaTwo) /
						sqrt(::xolotl::core::pi * twoSigmaTwo);

					// Update the rate
					f = this->_coefs(0, 0, 0, 0) * rate * gauss;

					// The large bubble concentration decreases
					if (this->_reactants[0] >= numClusters) {
						Kokkos::atomic_sub(
							&values(this->_connEntries[0][0][0][0]), f * stdC);
						Kokkos::atomic_sub(
							&values(this->_connEntries[0][0][1][0]), f * bC);
					}
					else {
						Kokkos::atomic_sub(
							&values(this->_connEntries[1][0][1][0]), f * stdC);
						Kokkos::atomic_sub(
							&values(this->_connEntries[1][0][0][0]), f * bC);
					}
					// The product concentration increases
					if (this->_reactants[0] >= numClusters) {
						Kokkos::atomic_add(
							&values(this->_connEntries[2][0][0][0]), f * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[2][0][1][0]), f * bC);
					}
					else {
						Kokkos::atomic_add(
							&values(this->_connEntries[2][0][1][0]), f * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[2][0][0][0]), f * bC);
					}

					// The V size decreases even more
					f = this->_coefs(0, 0, 0, 0) * rate *
						prComp[Species::V] *
						gauss;
					if (this->_reactants[0] >= numClusters) {
						Kokkos::atomic_sub(
							&values(this->_connEntries[0][1][0][0]), f * stdC);
						Kokkos::atomic_sub(
							&values(this->_connEntries[0][1][1][0]), f * bC);
					}
					else {
						Kokkos::atomic_sub(
							&values(this->_connEntries[1][1][1][0]), f * stdC);
						Kokkos::atomic_sub(
							&values(this->_connEntries[1][1][0][0]), f * bC);
					}
			}
				// In every case

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
						&values(this->_connEntries[0][1][0][0]), f * stdC);
					Kokkos::atomic_sub(
						&values(this->_connEntries[0][1][1][0]), f * bC);
				}
				else {
					Kokkos::atomic_sub(
						&values(this->_connEntries[1][1][1][0]), f * stdC);
					Kokkos::atomic_sub(
						&values(this->_connEntries[1][1][0][0]), f * bC);
				}
			}
*/
		// Helium case
		if (comp[Species::He] > 0) {
			// It should always be a bubble here

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
		
/*			// Compute the average concentrations
				auto conc = concentrations(ssbmId);
				auto avV = concentrations(ssbmId + 1) / conc;
				auto avHe = concentrations(ssbmId + 2) / conc;
				if (conc == 0.0) {
					avV = 0.0;
					avHe = 0.0;
				}

				// Sigmoid
				double x = ((comp[Species::He] + avHe) /
							   fe::getMaxHePerV(avV,avHe/avV)) -
					1.0;
				double sigmo = 1.0 / (1.0 + std::exp(-50.0 * x));
				double nI = 1.0;

				// The other product increases
				f = this->_coefs(0, 0, 0, 0) * rate * nI;
				if (this->_products[1] != Superclass::invalidIndex) {
					if (this->_reactants[0] >= numClusters) {
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][0][0]),
							f * sigmo * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][1][0]),
							f * sigmo * bC);
					}
					else {
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][1][0]),
							f * sigmo * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][0][0]),
							f * sigmo * bC);
					}
					// The average V increases
					Composition prodComp(
						this->_clusterData->getCluster(this->_products[1])
							.getRegion()
							.getOrigin());
					f = this->_coefs(0, 0, 0, 0) * rate * nI *
						prodComp[Species::I];
					if (this->_reactants[0] >= numClusters) {
						Kokkos::atomic_add(
							&values(this->_connEntries[0][2][0][0]),
							f * sigmo * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[0][2][1][0]),
							f * sigmo * bC);
					}
					else {
						Kokkos::atomic_add(
							&values(this->_connEntries[1][2][1][0]),
							f * sigmo * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[1][2][0][0]),
							f * sigmo * bC);
					}
				}
				else {
					f = this->_coefs(0, 0, 0, 0) * rate * nI;
					if (this->_reactants[0] >= numClusters) {
						Kokkos::atomic_add(
							&values(this->_connEntries[0][2][0][0]),
							f * sigmo * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[0][2][1][0]),
							f * sigmo * bC);
					}
					else {
						Kokkos::atomic_add(
							&values(this->_connEntries[1][2][1][0]),
							f * sigmo * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[1][2][0][0]),
							f * sigmo * bC);
					}
				}
*/
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

			// The large bubble increases, as well as average V
			Kokkos::atomic_add(
				&values(this->_connEntries[2][0][0][0]), f * cR2);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][0][1][0]), f * cR1);
			f = this->_coefs(0, 0, 0, 0) * rate *
				(comp1[Species::V] + comp2[Species::V]);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][1][0][0]), f * cR2);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][1][1][0]), f * cR1);
		}

		// Interstitial case
		if (orig1.isOnAxis(Species::I) or orig2.isOnAxis(Species::I)) {
			// Both reactants decrease
			Kokkos::atomic_sub(
				&values(this->_connEntries[0][0][0][0]), f * cR2);
			Kokkos::atomic_sub(
				&values(this->_connEntries[1][0][0][0]), f * cR2);
			Kokkos::atomic_sub(
				&values(this->_connEntries[0][0][1][0]), f * cR1);
			Kokkos::atomic_sub(
				&values(this->_connEntries[1][0][1][0]), f * cR1);

			// The large bubble increases, as well as average I
			Kokkos::atomic_add(
				&values(this->_connEntries[2][0][0][0]), f * cR2);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][0][1][0]), f * cR1);
			f = this->_coefs(0, 0, 0, 0) * rate *
				(comp1[Species::I] + comp2[Species::I]);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][1][0][0]), f * cR2);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][1][1][0]), f * cR1);
		}

		// Helium case
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

			// The large bubble increases, as well as average He
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
		}
	}
}


KOKKOS_INLINE_FUNCTION
double
FeDissociationReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return fe::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
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
			auto amtHe = comp[Species::He], amtV = comp[Species::V];
			double omega = this->_clusterData->atomicVolume();
			double T = this->_clusterData->temperature(0);
			constexpr double k_B = ::xolotl::core::kBoltzmann;
			double dg = 0.3135 * (0.8542 - 0.03996 * log(9.16 / T));
			double eta = ::xolotl::core::pi * dg * dg * dg * (double)amtHe /
				(6.0 * omega * (double)amtV);
			double z = (1.0 + eta + eta * eta * (1.0 - eta)) /
				((1.0 - eta) * (1.0 - eta) * (1.0 - eta));
			double p = (double)amtHe * z * k_B * T / (double)amtV;
			// HeV -> V
			if (prod1Comp.isOnAxis(Species::V) ||
				prod2Comp.isOnAxis(Species::V)) {
				return 5.0;
				if (amtV == 1 and amtHe < 9)
					be = bubbleV1Table[amtHe - 1];
				else if (amtV < 5 and amtHe < 5) {
					be = bubbleVTable[amtV - 2][amtHe - 1];
				}
				else
					be = 1.73 -
						2.59 *
							(pow((double)amtV, 2.0 / 3.0) -
								pow((double)amtV - 1.0, 2.0 / 3.0)) +
						2.50 * log( 1 + ((double)amtHe) /((double)amtV));
			}
			// HeV -> I
			if (prod1Comp.isOnAxis(Species::I) ||
				prod2Comp.isOnAxis(Species::I)) {
				if (amtV < 4 and amtHe < 5)
					be = bubbleITable[amtV - 1][amtHe - 1];
				else
					be = 4.88 +
						2.59 *
							(pow((double)amtV, 2.0 / 3.0) -
								pow((double)amtV - 1.0, 2.0 / 3.0)) -
						2.50 * log( 1 + ((double)amtHe) /((double)amtV));
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
		double omega = this->_clusterData->atomicVolume();
		double T = this->_clusterData->temperature(0);
		constexpr double k_B = ::xolotl::core::kBoltzmann;
		double dg = 0.3135 * (0.8542 - 0.03996 * log(9.16 / T));
		double eta = ::xolotl::core::pi * dg * dg * dg * (double)amtHe /
			(6.0 * omega * (double)amtV);
		double z = (1.0 + eta + eta * eta * (1.0 - eta)) /
			((1.0 - eta) * (1.0 - eta) * (1.0 - eta));
		double p = (double)amtHe * z * k_B * T / (double)amtV;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			return 5.0;
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
