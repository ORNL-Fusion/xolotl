#pragma once

#include <xolotl/core/network/PSIClusterGenerator.h>
#include <xolotl/core/network/impl/Reaction.tpp>
#include <xolotl/core/network/impl/SinkReaction.tpp>
#include <xolotl/core/network/impl/TrapMutationReaction.tpp>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace psi
{
template <typename TRegion>
KOKKOS_INLINE_FUNCTION
double
getRate(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg, const double r0,
	const double r1, const double dc0, const double dc1)
{
	constexpr double pi = ::xolotl::core::pi;

	double kPlus = 4.0 * pi * (r0 + r1) * (dc0 + dc1);

	return kPlus;
}
} // namespace psi

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSIProductionReaction<TSpeciesEnum>::getRateForProduction(IndexType gridIndex)
{
	// Standard case
	if (not isLargeBubbleReaction) {
		auto cl0 = this->_clusterData->getCluster(this->_reactants[0]);
		auto cl1 = this->_clusterData->getCluster(this->_reactants[1]);

		double r0 = cl0.getReactionRadius();
		double r1 = cl1.getReactionRadius();

		double dc0 = cl0.getDiffusionCoefficient(gridIndex);
		double dc1 = cl1.getDiffusionCoefficient(gridIndex);

		return psi::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
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

	return psi::getRate(cl0Reg, cl1Reg, r0, r1, dc0, dc1);
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSIDissociationReaction<TSpeciesEnum>::getRateForProduction(IndexType gridIndex)
{
	// Standard case
	if (not isLargeBubbleReaction) {
		auto cl0 = this->_clusterData->getCluster(this->_products[0]);
		auto cl1 = this->_clusterData->getCluster(this->_products[1]);

		double r0 = cl0.getReactionRadius();
		double r1 = cl1.getReactionRadius();

		double dc0 = cl0.getDiffusionCoefficient(gridIndex);
		double dc1 = cl1.getDiffusionCoefficient(gridIndex);

		return psi::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
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

	return psi::getRate(cl0Reg, cl1Reg, r0, r1, dc0, dc1);
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
void
PSIProductionReaction<TSpeciesEnum>::computeCoefficients()
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

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
void
PSIProductionReaction<TSpeciesEnum>::computeFlux(
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
			auto avV = concentrations(ssbmId + 2) / concentrations(ssbmId);
			auto avH = concentrations(ssbmId + 1) / concentrations(ssbmId);
			if (concentrations(ssbmId) == 0.0) {
				avV = 0.0;
				avH = 0.0;
			}

			// Get the largest cluster
			double largestVSize = this->_clusterData->maxVSize();

			// Threshold contribution
			double sigmo = util::computeSigmoid(avV, largestVSize, 2.0);

			// H/V contribution
			double maxH = psi::getMaxHPerV(util::max(avV - 1.0, 1.0),
				this->_clusterData->latticeParameter(),
				cl.getTemperature(gridIndex));
			sigmo *= 1.0 - util::computeSigmoid(avH, maxH, 2.0);

			// The standard cluster always loses the flux
			Kokkos::atomic_sub(&fluxes[stdClusterId], f * sigmo);

			// The V size decreases
			Kokkos::atomic_sub(
				&fluxes[ssbmId + 2], f * comp[Species::I] * sigmo);
		}

		// H case
		// H_k + B -> B
		if constexpr (psi::hasDeuterium<Species>) {
			if (comp[Species::D] > 0) {
				// The standard cluster always loses the flux
				Kokkos::atomic_sub(&fluxes[stdClusterId], f);

				// The average H increases
				Kokkos::atomic_add(&fluxes[ssbmId + 1], f * comp[Species::D]);

				// Trap mutation case
				// Compute the average concentrations
				auto conc = concentrations(ssbmId);
				auto avV = concentrations(ssbmId + 2) / conc;
				auto avH = concentrations(ssbmId + 1) / conc;
				if (conc == 0.0) {
					avV = 0.0;
					avH = 0.0;
				}

				// Sigmoid
				double maxH = psi::getMaxHPerV(avV,
					this->_clusterData->latticeParameter(),
					cl.getTemperature(gridIndex));
				double sigmo =
					util::computeSigmoid(comp[Species::D] + avH, maxH, 2.0);
				if (maxH == 0.0)
					sigmo = 1.0;
				if (avH < 1.0)
					sigmo = 0.0;
				if (avV <= 0.0 or conc < 1.0e-16)
					sigmo = 0.0;

				// The other product increases
				if (this->_products[1] != Superclass::invalidIndex) {
					Kokkos::atomic_add(&fluxes[this->_products[1]], f * sigmo);

					// The average V increases
					Composition prodComp(
						this->_clusterData->getCluster(this->_products[1])
							.getRegion()
							.getOrigin());
					Kokkos::atomic_add(
						&fluxes[ssbmId + 2], f * prodComp[Species::I] * sigmo);
				}
				else {
					// I is not explicitely modeled
					Kokkos::atomic_add(&fluxes[ssbmId + 2], f * sigmo);
				}
			}
		}

		// He case
		// He_k + B -> B
		// Structurally parallel to the H block above, retargeted to He. 
		if (comp[Species::He] > 0) {
			// The standard cluster always loses the flux
			Kokkos::atomic_sub(&fluxes[stdClusterId], f);

			// The average He increases
			Kokkos::atomic_add(&fluxes[ssbmId + 1], f * comp[Species::He]);

			// Trap mutation case
			// Compute the average concentrations. Guard before dividing:
			// Newton produces tiny and occasionally negative concentrations
			// on intermediate iterates, and dividing first raises FE flags
			// even when the result is later overwritten.
			auto conc = concentrations(ssbmId);
			double avV = 0.0;
			double avHe = 0.0;
			if (conc > 1.0e-16) {
				avV = concentrations(ssbmId + 2) / conc;
				avHe = concentrations(ssbmId + 1) / conc;
			}

			// Sigmoid
			double maxHe = psi::getMaxHePerVCont(avV,
				this->_clusterData->latticeParameter(),
				cl.getTemperature(gridIndex));
			double sigmo =
				util::computeSigmoid(comp[Species::He] + avHe, maxHe, 2.0);
			if (maxHe == 0.0)
				sigmo = 1.0;
			if (avHe < 1.0)
				sigmo = 0.0;
			if (avV <= 0.0 or conc < 1.0e-16)
				sigmo = 0.0;

			// The other product increases
			if (this->_products[1] != Superclass::invalidIndex) {
				Kokkos::atomic_add(&fluxes[this->_products[1]], f * sigmo);

				// The average V increases
				Composition prodComp(
					this->_clusterData->getCluster(this->_products[1])
						.getRegion()
						.getOrigin());
				Kokkos::atomic_add(
					&fluxes[ssbmId + 2], f * prodComp[Species::I] * sigmo);
			}
			else {
				// I is not currently explicitly modeled
				Kokkos::atomic_add(&fluxes[ssbmId + 2], f * sigmo);
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

		double f = this->_coefs(0, 0, 0, 0) * cR1 * cR2 * rate;

		// X_a + Y_b -> B
		//
		// SINGLE DISPATCH. The original code had one independent `if` per
		// species axis. For a reaction such as He_a + V_b -> B, both the V
		// test and the He test are true, so both blocks fired: the reactants
		// were subtracted twice and the product added twice. 
		//
		// Here the reaction is handled exactly once and the moments are
		// accumulated from the same totals.
		bool isBubbleForming = orig1.isOnAxis(Species::V) or
			orig2.isOnAxis(Species::V) or orig1.isOnAxis(Species::He) or
			orig2.isOnAxis(Species::He);
		if constexpr (psi::hasDeuterium<Species>) {
			isBubbleForming = isBubbleForming or
				orig1.isOnAxis(Species::D) or orig2.isOnAxis(Species::D);
		}

		if (isBubbleForming) {
			// Both reactants decrease
			Kokkos::atomic_sub(&fluxes[this->_reactants[0]], f);
			Kokkos::atomic_sub(&fluxes[this->_reactants[1]], f);

			// The large bubble increases
			Kokkos::atomic_add(&fluxes[this->_products[0]], f);

			// The average V increases
			auto totalVSize = comp1[Species::V] + comp2[Species::V];
			Kokkos::atomic_add(
				&fluxes[this->_products[0] + 2], f * totalVSize);

			// The average gas content increases. He and the H isotopes share
			// the single gas moment slot (see PSIReaction.h), so they are
			// summed here.
			auto totalGasSize = comp1[Species::He] + comp2[Species::He];
			if constexpr (psi::hasDeuterium<Species>) {
				totalGasSize += comp1[Species::D] + comp2[Species::D];
			}
			Kokkos::atomic_add(
				&fluxes[this->_products[0] + 1], f * totalGasSize);
		}
	}
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
void
PSIProductionReaction<TSpeciesEnum>::computePartialDerivatives(
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
	//
	// _connEntries is declared
	//     util::Array<IndexType, 4, 1 + nMomentIds, 2, 1 + nMomentIds>
	// so both the second and the fourth index range over 1 + nMomentIds with
	// the same layout: 0 is the concentration, 1..nMomentIds are the moments.
	// The gas (He) moment is always slot 1; the V moment is always the last
	// slot, nMomentIds.
	//
	// For the full He/D/V network nMomentIds == 3, so V is slot 3 and the
	// hardcoded literals that used to appear here happened to be right. For a
	// pure-He network nMomentIds == 2, so V is slot 2 and a literal 3 reads
	// past the end of the array.
	constexpr auto gasMomId = 1;
	constexpr auto vMomId = Superclass::nMomentIds;

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
					&values(this->_connEntries[0][vMomId][0][0]), f * stdC);
				Kokkos::atomic_add(
					&values(this->_connEntries[0][vMomId][1][0]), f * bC);
			}
			else {
				Kokkos::atomic_add(
					&values(this->_connEntries[1][vMomId][1][0]), f * stdC);
				Kokkos::atomic_add(
					&values(this->_connEntries[1][vMomId][0][0]), f * bC);
			}
		}
		// Interstitial case
		if (comp[Species::I] > 0) {
			// I_k + B -> B
			// The reaction can happen only if the average V is larger than the
			// threshold and if the H/V ratio is realistic
			auto avV = concentrations(ssbmId + 2) / concentrations(ssbmId);
			auto avH = concentrations(ssbmId + 1) / concentrations(ssbmId);
			if (concentrations(ssbmId) == 0.0) {
				avV = 0.0;
				avH = 0.0;
			}

			// Get the largest cluster
			double largestVSize = this->_clusterData->maxVSize();

			// Threshold contribution
			double sigmo = util::computeSigmoid(avV, largestVSize, 2.0);

			// H/V contribution
			double maxH = psi::getMaxHPerV(util::max(avV - 1.0, 1.0),
				this->_clusterData->latticeParameter(),
				cl.getTemperature(gridIndex));
			sigmo *= 1.0 - util::computeSigmoid(avH, maxH, 2.0);

			// The standard cluster always loses the flux
			f = this->_coefs(0, 0, 0, 0) * rate * sigmo;
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
			f = this->_coefs(0, 0, 0, 0) * rate * comp[Species::I] * sigmo;
			if (this->_reactants[0] >= numClusters) {
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][vMomId][0][0]), f * stdC);
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][vMomId][1][0]), f * bC);
			}
			else {
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][vMomId][1][0]), f * stdC);
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][vMomId][0][0]), f * bC);
			}
		}

		// H case
		if constexpr (psi::hasDeuterium<Species>) {
			if (comp[Species::D] > 0) {
				// H_k + B -> B

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

				// The H size increases
				f = this->_coefs(0, 0, 0, 0) * rate * comp[Species::D];
				if (this->_reactants[0] >= numClusters) {
					Kokkos::atomic_add(
						&values(this->_connEntries[0][gasMomId][0][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[0][gasMomId][1][0]), f * bC);
				}
				else {
					Kokkos::atomic_add(
						&values(this->_connEntries[1][gasMomId][1][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[1][gasMomId][0][0]), f * bC);
				}

				// Compute the average concentrations
				auto conc = concentrations(ssbmId);
				auto avV = concentrations(ssbmId + 2) / conc;
				auto avH = concentrations(ssbmId + 1) / conc;
				if (conc == 0.0) {
					avV = 0.0;
					avH = 0.0;
				}

				// Sigmoid for trap mutation
				// H_k + B -> B + I
				double maxH = psi::getMaxHPerV(avV,
					this->_clusterData->latticeParameter(),
					cl.getTemperature(gridIndex));
				double sigmo =
					util::computeSigmoid(comp[Species::D] + avH, maxH, 2.0);
				if (maxH == 0.0)
					sigmo = 1.0;
				if (avH < 1.0)
					sigmo = 0.0;
				if (avV <= 0.0 or bC < 1.0e-16)
					sigmo = 0.0;

				// The other product increases
				f = this->_coefs(0, 0, 0, 0) * rate * sigmo;
				if (this->_products[1] != Superclass::invalidIndex) {
					// Compute the terms depending on average H and V
					double w = psi::getMaxHPerV(avV,
						this->_clusterData->latticeParameter(),
						cl.getTemperature(gridIndex));
					double wdV = psi::getMaxHPerVdV(avV,
						this->_clusterData->latticeParameter(),
						cl.getTemperature(gridIndex));
					double y = -2.0 * (comp[Species::D] + avH) + 2.0 * w;
					double ydH = -2.0;
					double ydV = 2.0 * wdV;
					double z = 1.0 + exp(y);
					double zdH = ydH * exp(y);
					double zdV = ydV * exp(y);
					double g = 1.0 / z;
					double gdH = -zdH / (z * z);
					double gdV = -zdV / (z * z);
					if (g == 0.0) {
						gdH = 0.0;
						gdV = 0.0;
					}
					double fdH = this->_coefs(0, 0, 0, 0) * rate * gdH * stdC;
					double fdV = this->_coefs(0, 0, 0, 0) * rate * gdV * stdC;
					if (avV <= 0.0 or bC < 1.0e-16) {
						fdH = 0.0;
						fdV = 0.0;
					}
					if (this->_reactants[0] >= numClusters) {
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][0][0]), f * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][0][gasMomId]), fdH);
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][0][vMomId]), fdV);
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][1][0]), f * bC);
					}
					else {
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][1][0]), f * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][1][gasMomId]), fdH);
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][1][vMomId]), fdV);
						Kokkos::atomic_add(
							&values(this->_connEntries[3][0][0][0]), f * bC);
					}
					// The average V increases
					Composition prodComp(
						this->_clusterData->getCluster(this->_products[1])
							.getRegion()
							.getOrigin());
					f = this->_coefs(0, 0, 0, 0) * rate * prodComp[Species::I] *
						sigmo;
					if (this->_reactants[0] >= numClusters) {
						Kokkos::atomic_add(
							&values(this->_connEntries[0][vMomId][0][0]), f * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[0][vMomId][0][gasMomId]), fdH);
						Kokkos::atomic_add(
							&values(this->_connEntries[0][vMomId][0][vMomId]), fdV);
						Kokkos::atomic_add(
							&values(this->_connEntries[0][vMomId][1][0]), f * bC);
					}
					else {
						Kokkos::atomic_add(
							&values(this->_connEntries[1][vMomId][1][0]), f * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[1][vMomId][1][gasMomId]), fdH);
						Kokkos::atomic_add(
							&values(this->_connEntries[1][vMomId][1][vMomId]), fdV);
						Kokkos::atomic_add(
							&values(this->_connEntries[1][vMomId][0][0]), f * bC);
					}
				}
				else {
					// The average V increases but I is not explicitly modeled
					// here
					f = this->_coefs(0, 0, 0, 0) * rate * sigmo;
					if (this->_reactants[0] >= numClusters) {
						Kokkos::atomic_add(
							&values(this->_connEntries[0][vMomId][0][0]), f * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[0][vMomId][1][0]), f * bC);
					}
					else {
						Kokkos::atomic_add(
							&values(this->_connEntries[1][vMomId][1][0]), f * stdC);
						Kokkos::atomic_add(
							&values(this->_connEntries[1][vMomId][0][0]), f * bC);
					}
				}
			}
		}

		// He case
		if (comp[Species::He] > 0) {
			// He_k + B -> B

			// Reset f: it is reassigned by each preceding species block, and
			// relying on the value left behind by another block would silently
			// carry that block's sigmoid into this one.
			f = this->_coefs(0, 0, 0, 0) * rate;

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
					&values(this->_connEntries[0][gasMomId][0][0]), f * stdC);
				Kokkos::atomic_add(
					&values(this->_connEntries[0][gasMomId][1][0]), f * bC);
			}
			else {
				Kokkos::atomic_add(
					&values(this->_connEntries[1][gasMomId][1][0]), f * stdC);
				Kokkos::atomic_add(
					&values(this->_connEntries[1][gasMomId][0][0]), f * bC);
			}

			// Compute the average concentrations (guard before dividing)
			double avV = 0.0;
			double avHe = 0.0;
			if (bC > 1.0e-16) {
				avV = concentrations(ssbmId + 2) / bC;
				avHe = concentrations(ssbmId + 1) / bC;
			}

			// Sigmoid for trap mutation
			// He_k + B -> B + I
			// These guards must match computeFlux. If the residual and
			// the Jacobian apply different guards they describe different
			// functions and Newton will not converge.
			double maxHe = psi::getMaxHePerVCont(avV,
				this->_clusterData->latticeParameter(),
				cl.getTemperature(gridIndex));
			double sigmo =
				util::computeSigmoid(comp[Species::He] + avHe, maxHe, 2.0);
			if (maxHe == 0.0)
				sigmo = 1.0;
			if (avHe < 1.0)
				sigmo = 0.0;
			if (avV <= 0.0 or bC < 1.0e-16)
				sigmo = 0.0;

			// The other product increases
			f = this->_coefs(0, 0, 0, 0) * rate * sigmo;
			if (this->_products[1] != Superclass::invalidIndex) {
				// Chain rule through the sigmoid. The flux depends on avHe and
				// avV, which are moment/concentration, so there are derivative
				// paths with respect to both moment DOFs. Omitting these
				// leaves the Jacobian inconsistent with the residual.
				//   g(y) = 1/(1+exp(y)),
				//   y = -sharp*(He_k + <He>) + sharp*w,  w = maxHe(<V>)
				double w = psi::getMaxHePerVCont(avV,
					this->_clusterData->latticeParameter(),
					cl.getTemperature(gridIndex));
				double wdV = psi::getMaxHePerVdV(avV,
					this->_clusterData->latticeParameter(),
					cl.getTemperature(gridIndex));
				double y = -2.0 * (comp[Species::He] + avHe) + 2.0 * w;
				double ydHe = -2.0;
				double ydV = 2.0 * wdV;
				double z = 1.0 + exp(y);
				double zdHe = ydHe * exp(y);
				double zdV = ydV * exp(y);
				double g = 1.0 / z;
				double gdHe = -zdHe / (z * z);
				double gdV = -zdV / (z * z);
				if (g == 0.0) {
					gdHe = 0.0;
					gdV = 0.0;
				}
				double fdHe = this->_coefs(0, 0, 0, 0) * rate * gdHe * stdC;
				double fdV = this->_coefs(0, 0, 0, 0) * rate * gdV * stdC;
				if (avV <= 0.0 or bC < 1.0e-16) {
					fdHe = 0.0;
					fdV = 0.0;
				}
				if (this->_reactants[0] >= numClusters) {
					Kokkos::atomic_add(
						&values(this->_connEntries[3][0][0][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[3][0][0][gasMomId]), fdHe);
					Kokkos::atomic_add(
						&values(this->_connEntries[3][0][0][vMomId]), fdV);
					Kokkos::atomic_add(
						&values(this->_connEntries[3][0][1][0]), f * bC);
				}
				else {
					Kokkos::atomic_add(
						&values(this->_connEntries[3][0][1][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[3][0][1][gasMomId]), fdHe);
					Kokkos::atomic_add(
						&values(this->_connEntries[3][0][1][vMomId]), fdV);
					Kokkos::atomic_add(
						&values(this->_connEntries[3][0][0][0]), f * bC);
				}
				// The average V increases
				Composition prodComp(
					this->_clusterData->getCluster(this->_products[1])
						.getRegion()
						.getOrigin());
				f = this->_coefs(0, 0, 0, 0) * rate * prodComp[Species::I] *
					sigmo;
				if (this->_reactants[0] >= numClusters) {
					Kokkos::atomic_add(
						&values(this->_connEntries[0][vMomId][0][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[0][vMomId][0][gasMomId]),
						fdHe);
					Kokkos::atomic_add(
						&values(this->_connEntries[0][vMomId][0][vMomId]), fdV);
					Kokkos::atomic_add(
						&values(this->_connEntries[0][vMomId][1][0]), f * bC);
				}
				else {
					Kokkos::atomic_add(
						&values(this->_connEntries[1][vMomId][1][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[1][vMomId][1][gasMomId]),
						fdHe);
					Kokkos::atomic_add(
						&values(this->_connEntries[1][vMomId][1][vMomId]), fdV);
					Kokkos::atomic_add(
						&values(this->_connEntries[1][vMomId][0][0]), f * bC);
				}
			}
			else {
				// The average V increases but I is not explicitly modeled
				f = this->_coefs(0, 0, 0, 0) * rate * sigmo;
				if (this->_reactants[0] >= numClusters) {
					Kokkos::atomic_add(
						&values(this->_connEntries[0][vMomId][0][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[0][vMomId][1][0]), f * bC);
				}
				else {
					Kokkos::atomic_add(
						&values(this->_connEntries[1][vMomId][1][0]), f * stdC);
					Kokkos::atomic_add(
						&values(this->_connEntries[1][vMomId][0][0]), f * bC);
				}
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

		// X_a + Y_b -> B
		// Single dispatch, mirroring computeFlux. See the comment there: one
		// independent `if` per species axis double-counted any reaction whose
		// reactants sat on two of the tested axes.
		bool isBubbleForming = orig1.isOnAxis(Species::V) or
			orig2.isOnAxis(Species::V) or orig1.isOnAxis(Species::He) or
			orig2.isOnAxis(Species::He);
		if constexpr (psi::hasDeuterium<Species>) {
			isBubbleForming = isBubbleForming or
				orig1.isOnAxis(Species::D) or orig2.isOnAxis(Species::D);
		}

		if (isBubbleForming) {
			// Both reactants decrease
			Kokkos::atomic_sub(
				&values(this->_connEntries[0][0][0][0]), f * cR2);
			Kokkos::atomic_sub(
				&values(this->_connEntries[1][0][0][0]), f * cR2);
			Kokkos::atomic_sub(
				&values(this->_connEntries[0][0][1][0]), f * cR1);
			Kokkos::atomic_sub(
				&values(this->_connEntries[1][0][1][0]), f * cR1);

			// The large bubble increases
			Kokkos::atomic_add(
				&values(this->_connEntries[2][0][0][0]), f * cR2);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][0][1][0]), f * cR1);

			// The average gas content increases (He and H isotopes share the
			// single gas moment slot)
			auto totalGasSize = comp1[Species::He] + comp2[Species::He];
			if constexpr (psi::hasDeuterium<Species>) {
				totalGasSize += comp1[Species::D] + comp2[Species::D];
			}
			f = this->_coefs(0, 0, 0, 0) * rate * totalGasSize;
			Kokkos::atomic_add(
				&values(this->_connEntries[2][gasMomId][0][0]), f * cR2);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][gasMomId][1][0]), f * cR1);

			// The average V increases
			f = this->_coefs(0, 0, 0, 0) * rate *
				(comp1[Species::V] + comp2[Species::V]);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][vMomId][0][0]), f * cR2);
			Kokkos::atomic_add(
				&values(this->_connEntries[2][vMomId][1][0]), f * cR1);
		}
	}
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
void
PSIDissociationReaction<TSpeciesEnum>::computeCoefficients()
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

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
void
PSIDissociationReaction<TSpeciesEnum>::computeFlux(
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
	double E_b = 2.15; // H in Void

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
	auto ssbmId = (this->_products[0] >= numClusters) ? this->_products[0] :
														this->_products[1];

	// Compute the flux
	double f = this->_coefs(0, 0, 0, 0) * rate * concentrations(stdClusterId) *
		concentrations(ssbmId);

	// Large bubble is the reactant
	if (this->_reactant >= numClusters) {
		// H case
		// B -> H_1 + B
		if constexpr (psi::hasDeuterium<Species>) {
			if (comp[Species::D] > 0) {
				// The standard cluster always gains the flux
				Kokkos::atomic_add(&fluxes[stdClusterId], f);

				// The average H decrease
				Kokkos::atomic_sub(&fluxes[ssbmId + 1], f * comp[Species::D]);
			}
		}
	}
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
void
PSIDissociationReaction<TSpeciesEnum>::computePartialDerivatives(
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
	double E_b = 2.15; // H in Void

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
	auto ssbmId = (this->_products[0] >= numClusters) ? this->_products[0] :
														this->_products[1];

	// Get the concentration
	auto bC = concentrations(ssbmId);

	// Large bubble is the reactant
	if (this->_reactant >= numClusters) {
		// H case
		// B -> H_1 + B
		if constexpr (psi::hasDeuterium<Species>) {
			if (comp[Species::D] > 0) {
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
				f = this->_coefs(0, 0, 0, 0) * rate * comp[Species::D];
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
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSIDissociationReaction<TSpeciesEnum>::computeBindingEnergy(double time)
{
	using psi::hasDeuterium;
	using psi::hasTritium;

	using NetworkType = typename Superclass::NetworkType;

	constexpr double beTableV1[10][7] = {
		// H:  1     2     3     4     5     6      // He:
		{0.0, 1.38, 1.34, 1.08, 1.04, 0.92, 0.71}, // 0
		{0.0, 1.00, 0.95, 0.90, 0.88, 0.80, 0.60}, // 1
		{0.0, 0.96, 0.92, 0.85, 0.84, 0.83, 0.50}, // 2
		{0.0, 0.86, 0.81, 0.69, 0.64, 0.65, 0.50}, // 3
		{0.0, 0.83, 0.80, 0.65, 0.60, 0.60, 0.55}, // 4
		{0.0, 0.83, 0.80, 0.60, 0.50, 0.50, 0.50}, // 5
		{0.0, 0.80, 0.70, 0.60, 0.50, 0.50, 0.50}, // 6
		{0.0, 0.80, 0.75, 0.65, 0.55, 0.55, 0.45}, // 7
		{0.0, 0.80, 0.80, 0.70, 0.65, 0.60, 0.55}, // 8
		{0.0, 0.80, 0.80, 0.75, 0.70, 0.65, 0.60}, // 9
	};

	constexpr double beTableV2[15][12] = {
		// H:  1     2     3     4     5     6     7     8     9     10    11 //
		// He:
		{0.0, 1.60, 1.60, 1.39, 1.30, 1.15, 1.15, 1.05, 1.05, 0.85, 0.80,
			0.70}, // 0
		{0.0, 1.42, 1.45, 1.35, 1.04, 1.00, 1.00, 1.00, 1.00, 0.8, 0.75,
			0.70}, // 1
		{0.0, 1.18, 1.15, 1.15, 1.21, 1.00, 1.00, 1.00, 0.95, 0.80, 0.75,
			0.70}, // 2
		{0.0, 1.25, 1.16, 1.04, 1.05, 0.95, 1.00, 0.95, 0.95, 0.75, 0.75,
			0.70}, // 3
		{0.0, 1.15, 1.11, 1.04, 1.04, 0.95, 1.00, 0.95, 0.95, 0.75, 0.70,
			0.70}, // 4
		{0.0, 1.15, 1.10, 1.04, 1.04, 0.95, 0.95, 0.95, 0.95, 0.75, 0.70,
			0.70}, // 5
		{0.0, 1.15, 1.10, 1.04, 1.04, 0.95, 0.95, 0.95, 0.90, 0.75, 0.70,
			0.65}, // 6
		{0.0, 1.10, 1.05, 1.00, 1.00, 0.95, 0.95, 0.95, 0.90, 0.70, 0.70,
			0.65}, // 7
		{0.0, 1.10, 1.05, 1.00, 1.00, 0.95, 0.95, 0.90, 0.90, 0.70, 0.70,
			0.65}, // 8
		{0.0, 1.10, 1.05, 1.00, 1.00, 0.90, 0.90, 0.90, 0.90, 0.70, 0.70,
			0.65}, // 9
		{0.0, 1.05, 1.00, 0.95, 0.95, 0.90, 0.90, 0.90, 0.85, 0.70, 0.65,
			0.65}, // 10
		{0.0, 1.00, 1.00, 0.95, 0.95, 0.90, 0.90, 0.90, 0.85, 0.70, 0.65,
			0.65}, // 11
		{0.0, 1.00, 0.95, 0.95, 0.90, 0.90, 0.90, 0.85, 0.85, 0.65, 0.65,
			0.60}, // 12
		{0.0, 0.95, 0.95, 0.90, 0.90, 0.90, 0.90, 0.85, 0.85, 0.65, 0.65,
			0.60}, // 13
		{0.0, 0.95, 0.95, 0.90, 0.90, 0.85, 0.85, 0.85, 0.75, 0.65, 0.65,
			0.60}, // 14
	};

	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;

	double be = 0.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);
	auto prod1 = this->_clusterData->getCluster(this->_products[0]);
	auto prod2 = this->_clusterData->getCluster(this->_products[1]);

	auto clReg = cl.getRegion();
	auto prod1Reg = prod1.getRegion();
	auto prod2Reg = prod2.getRegion();
	bool useTable = false;
	if (clReg.isSimplex()) {
		if (prod1Reg.isSimplex()) {
			auto orig1 = prod1Reg.getOrigin();
			if constexpr (hasDeuterium<Species> && hasTritium<Species>) {
				if (orig1.isOnAxis(Species::D) || orig1.isOnAxis(Species::T)) {
					useTable = true;
				}
			}
		}
		if (prod2Reg.isSimplex()) {
			auto orig2 = prod2Reg.getOrigin();
			if constexpr (hasDeuterium<Species> && hasTritium<Species>) {
				if (orig2.isOnAxis(Species::D) || orig2.isOnAxis(Species::T)) {
					useTable = true;
				}
			}
		}
	}

	if constexpr (hasDeuterium<Species> && hasTritium<Species>) {
		if (useTable) {
			Composition comp(clReg.getOrigin());
			auto hAmount = comp[Species::D] + comp[Species::T];
			if (comp[Species::V] == 1) {
				be = beTableV1[comp[Species::He]][hAmount];
			}
			else if (comp[Species::V] == 2) {
				be = beTableV2[comp[Species::He]][hAmount];
			}
		}
	}

	if (be == 0.0) {
		// Special case for V
		auto orig1 = prod1Reg.getOrigin();
		auto orig2 = prod2Reg.getOrigin();
		Composition comp(clReg.getOrigin());
		AmountType lowerV = 16, higherV = 31;
		AmountType minV = 1;
		for (auto i = 1; i < higherV; i++) {
			auto maxHe = psi::getMaxHePerV(i, 4.0);
			if (comp[Species::He] > maxHe)
				minV = i;
		}
		lowerV = util::max(lowerV, minV + 2);
		if ((orig1.isOnAxis(Species::V) || orig2.isOnAxis(Species::V)) &&
			(comp[Species::V] >= lowerV && comp[Species::V] <= higherV)) {
			// Get the be at 16 and 30
			Composition HeVComp(clReg.getOrigin());
			HeVComp[Species::V] = lowerV;
			auto fe1 = PSIClusterGenerator<TSpeciesEnum>::getHeVFormationEnergy(
				HeVComp);
			HeVComp[Species::V] = lowerV - 1;
			auto fe2 = PSIClusterGenerator<TSpeciesEnum>::getHeVFormationEnergy(
				HeVComp);
			Composition vComp{};
			vComp[Species::V] = 1;
			auto fe3 =
				PSIClusterGenerator<TSpeciesEnum>::getHeVFormationEnergy(vComp);
			auto be1 = fe2 + fe3 - fe1;
			HeVComp[Species::V] = higherV;
			fe1 = PSIClusterGenerator<TSpeciesEnum>::getHeVFormationEnergy(
				HeVComp);
			HeVComp[Species::V] = higherV - 1;
			fe2 = PSIClusterGenerator<TSpeciesEnum>::getHeVFormationEnergy(
				HeVComp);
			auto be2 = fe2 + fe3 - fe1;
			if (higherV - lowerV < 4)
				be = be2;
			else
				be = be1 +
					(comp[Species::V] - lowerV) * (be2 - be1) /
						(higherV - lowerV);
		}
		else {
			be = prod1.getFormationEnergy() + prod2.getFormationEnergy() -
				cl.getFormationEnergy();
		}
	}

	return util::max(be, -5.0);
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSISinkReaction<TSpeciesEnum>::getSinkBias()
{
	return 1.0;
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSISinkReaction<TSpeciesEnum>::getSinkStrength()
{
	return this->_clusterData->sinkStrength();
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSITransformReaction<TSpeciesEnum>::getAppliedRate(IndexType gridIndex) const
{
	using NetworkType = typename Superclass::NetworkType;
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;

	// Get the radius of the cluster
	double radius = 0.0;
	if (isLargeBubbleReaction) {
		radius = this->_clusterData->bubbleAvRad();
	}
	else {
		auto cl = this->_clusterData->getCluster(this->_reactant);
		radius = cl.getReactionRadius();
	}

	// Get the current depth
	auto depth = this->_clusterData->getDepth();
	auto f = this->_clusterData->getFBursting();

	// TODO: change rate
	return f * util::computeSigmoid(radius, depth, 2.0);
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
void
PSITransformReaction<TSpeciesEnum>::computeCoefficients()
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

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
void
PSITransformReaction<TSpeciesEnum>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	// Standard case
	if (not isLargeBubbleReaction) {
		return Superclass::computeFlux(concentrations, fluxes, gridIndex);
	}

	auto rate = this->getAppliedRate(gridIndex);

	auto ssbmId = this->_product;
	auto avV = concentrations(ssbmId + 2) / concentrations(ssbmId);
	auto avH = concentrations(ssbmId + 1) / concentrations(ssbmId);
	if (concentrations(ssbmId) == 0.0) {
		avV = 0.0;
		avH = 0.0;
	}

	// Get the largest cluster
	double largestHSize = this->_clusterData->maxHSize();
	// Threshold contribution
	double sigmo = util::computeSigmoid(avH, largestHSize, 2.0);

	Kokkos::atomic_sub(
		&fluxes[ssbmId + 1], rate * concentrations[ssbmId + 1] * sigmo);
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
void
PSITransformReaction<TSpeciesEnum>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	// Standard case
	if (not isLargeBubbleReaction) {
		return Superclass::computePartialDerivatives(
			concentrations, values, gridIndex);
	}

	auto rate = this->getAppliedRate(gridIndex);

	auto ssbmId = this->_product;
	auto avV = concentrations(ssbmId + 2) / concentrations(ssbmId);
	auto avH = concentrations(ssbmId + 1) / concentrations(ssbmId);
	if (concentrations(ssbmId) == 0.0) {
		avV = 0.0;
		avH = 0.0;
	}

	// Get the largest cluster
	double largestHSize = this->_clusterData->maxHSize();
	// Threshold contribution
	double sigmo = util::computeSigmoid(avH, largestHSize, 2.0);

	Kokkos::atomic_sub(&values(this->_connEntries[1][2][0][0]), rate * sigmo);
}
} // namespace network
} // namespace core
} // namespace xolotl
