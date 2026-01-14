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
	auto cl0 = this->_clusterData->getCluster(this->_products[0]);
	auto cl1 = this->_clusterData->getCluster(this->_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return psi::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
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
		if (comp[Species::V] > 0) {
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
				// Only if the large void has a specific size
				auto avVoid =
					concentrations(ssbmId + 1) / concentrations(ssbmId);
				if (concentrations(ssbmId) == 0.0)
					avVoid = 0.0;

				// Get the product composition
				auto pr = this->_clusterData->getCluster(this->_products[0]);
				auto prReg = pr.getRegion();
				Composition prComp(prReg.getOrigin());

				// Target value for the reaction to happen
				double target = prComp[Species::V] + comp[Species::I];

				// Gaussian function around it
				double twoSigmaTwo = 0.1;
				gauss =
					exp(-(avVoid - target) * (avVoid - target) / twoSigmaTwo) /
					sqrt(::xolotl::core::pi * twoSigmaTwo);

				// The large bubble concentration decreases
				Kokkos::atomic_sub(&fluxes[ssbmId], f * gauss);
				// The product concentration increases
				Kokkos::atomic_add(&fluxes[this->_products[0]], f * gauss);

				// The V size decreases even more
				Kokkos::atomic_sub(
					&fluxes[ssbmId + 1], f * prComp[Species::V] * gauss);
			}

			// In every case

			// The standard cluster always loses the flux
			Kokkos::atomic_sub(&fluxes[stdClusterId], f * gauss);

			// The V size decreases
			Kokkos::atomic_sub(
				&fluxes[ssbmId + 1], f * comp[Species::I] * gauss);
		}

		// H case
		if constexpr (psi::hasDeuterium<Species>) {
			if (comp[Species::D] > 0) {
				// The standard cluster always loses the flux
				Kokkos::atomic_sub(&fluxes[stdClusterId], f);

				// The average H increases
				Kokkos::atomic_add(&fluxes[ssbmId + 2], f * comp[Species::D]);

				// Trap mutation case
				// Compute the average concentrations
				auto conc = concentrations(ssbmId);
				auto avV = concentrations(ssbmId + 1) / conc;
				auto avH = concentrations(ssbmId + 2) / conc;
				if (conc == 0.0) {
					avV = 0.0;
					avH = 0.0;
				}

				// Sigmoid
				double x = ((comp[Species::D] + avH) /
							   psi::getMaxHPerV(avV,
								   this->_clusterData->latticeParameter(),
								   cl.getTemperature(gridIndex))) -
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

		// Vacancy case
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

		// H case
		if constexpr (psi::hasDeuterium<Species>) {
			if (orig1.isOnAxis(Species::D) or orig2.isOnAxis(Species::D)) {
				// Compute the total size
				auto totalSize = comp1[Species::D] + comp2[Species::D];
				// Both reactants decrease
				Kokkos::atomic_sub(&fluxes[this->_reactants[0]], f);
				Kokkos::atomic_sub(&fluxes[this->_reactants[1]], f);
				// The large bubble increases, as well as average H
				Kokkos::atomic_add(&fluxes[this->_products[0]], f);
				Kokkos::atomic_add(
					&fluxes[this->_products[0] + 1], f * totalSize);
			}
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
		if (comp[Species::I] > 0) {
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
				auto pr = this->_clusterData->getCluster(this->_products[0]);
				auto prReg = pr.getRegion();
				Composition prComp(prReg.getOrigin());

				// Target value for the reaction to happen
				double target = prComp[Species::V] + comp[Species::I];

				// Gaussian function around it
				double twoSigmaTwo = 0.1;
				double gauss =
					exp(-(avVoid - target) * (avVoid - target) / twoSigmaTwo) /
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
				f = this->_coefs(0, 0, 0, 0) * rate * prComp[Species::V] *
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

		// H case
		if constexpr (psi::hasDeuterium<Species>) {
			if (comp[Species::D] > 0) {
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

				// The H size increases
				f = this->_coefs(0, 0, 0, 0) * rate * comp[Species::D];
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

				// Compute the average concentrations
				auto conc = concentrations(ssbmId);
				auto avV = concentrations(ssbmId + 1) / conc;
				auto avH = concentrations(ssbmId + 2) / conc;
				if (conc == 0.0) {
					avV = 0.0;
					avH = 0.0;
				}

				// Sigmoid
				double x = ((comp[Species::D] + avH) /
							   psi::getMaxHPerV(avV,
								   this->_clusterData->latticeParameter(),
								   cl.getTemperature(gridIndex))) -
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

		// H case
		if constexpr (psi::hasDeuterium<Species>) {
			if (orig1.isOnAxis(Species::D) or orig2.isOnAxis(Species::D)) {
				// Both reactants decrease
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][0][0][0]), f * cR2);
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][0][0][0]), f * cR2);
				Kokkos::atomic_sub(
					&values(this->_connEntries[0][0][1][0]), f * cR1);
				Kokkos::atomic_sub(
					&values(this->_connEntries[1][0][1][0]), f * cR1);

				// The large bubble increases, as well as average H
				Kokkos::atomic_add(
					&values(this->_connEntries[2][0][0][0]), f * cR2);
				Kokkos::atomic_add(
					&values(this->_connEntries[2][0][1][0]), f * cR1);
				f = this->_coefs(0, 0, 0, 0) * rate *
					(comp1[Species::D] + comp2[Species::D]);
				Kokkos::atomic_add(
					&values(this->_connEntries[2][2][0][0]), f * cR2);
				Kokkos::atomic_add(
					&values(this->_connEntries[2][2][1][0]), f * cR1);
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
} // namespace network
} // namespace core
} // namespace xolotl
