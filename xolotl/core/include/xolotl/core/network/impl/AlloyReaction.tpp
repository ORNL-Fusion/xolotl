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

	bool cl0IsSphere = (lo0[Species::He] > 0 or lo0[Species::V] > 0 or
			 lo0.isOnAxis(Species::I)),
		 cl1IsSphere = (lo1[Species::He] > 0 or lo1[Species::V] > 0 or
			 lo1.isOnAxis(Species::I));

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
		double temp = r1 / (3 * (r0 + rd));
		double alpha = 1.0 / (1.0 + temp * temp);
		double rateSpherical = 4.0 * pi * (r0 + r1 + rd);
		double rateToroidal =
			(4.0 * pi * pi * r1) / log(1 + (8 * r1) / (r0 + rd));

		// Calculate the capture efficiency (assuming only prismatic loops)
		if (cl1IsV)
			Pl = 0.78 / (p * p) + 0.66 * p - 0.44;
		else
			Pl = 0.70 / (p * p) + 0.78 * p - 0.47;

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

	// Create an array with all possible dislocation capture radii
	// rdCl = {(rdI for cl0, rdV for cl0), (rdI for cl1, rdV for cl1)}
	double rdCl[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

	// Standard case
	if (not isLargeBubbleReaction) {
		double r0 = cl0.getReactionRadius();
		double r1 = cl1.getReactionRadius();

		double dc0 = cl0.getDiffusionCoefficient(gridIndex);
		double dc1 = cl1.getDiffusionCoefficient(gridIndex);
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

	// Large bubble case
	const auto dummyRegion = Region(Composition{});

	double r0 = 0.0, r1 = 0.0, dc0 = 0.0, dc1 = 0.0;
	Region cl0Reg = dummyRegion, cl1Reg = dummyRegion;
	auto numClusters = this->_clusterData->numClusters;

	if (this->_reactants[0] >= numClusters) {
		auto i = (this->_reactants[0] - numClusters) / 2;
		switch (i) {
		// Bubble
		case 0:
			r0 = this->_clusterData->bubbleAvRad();
			cl0Reg[Species::V] = {1, 2};
			break;
		// Perfect V
		case 1:
			r0 = this->_clusterData->perfVAvRad();
			cl0Reg[Species::PerfectV] = {1, 2};
			rdCl[0][0] = 0.79;
			rdCl[0][1] = 1.59;
			break;
		// Faulted V
		case 2:
			r0 = this->_clusterData->faulVAvRad();
			cl0Reg[Species::FaultedV] = {1, 2};
			rdCl[0][0] = 0.79;
			rdCl[0][1] = 1.59;
			break;
		// Perfect I
		case 3:
			r0 = this->_clusterData->perfIAvRad();
			cl0Reg[Species::PerfectI] = {1, 2};
			rdCl[0][0] = 1.85;
			rdCl[0][1] = 1.04;
			break;
		// Faulted I
		case 4:
			r0 = this->_clusterData->faulIAvRad();
			cl0Reg[Species::FaultedI] = {1, 2};
			rdCl[0][0] = 1.85;
			rdCl[0][1] = 1.04;
			break;
		}
	}
	else {
		auto cl0 = this->_clusterData->getCluster(this->_reactants[0]);
		r0 = cl0.getReactionRadius();
		dc0 = cl0.getDiffusionCoefficient(gridIndex);
		cl0Reg = cl0.getRegion();
		rdCl[0][0] = this->_clusterData->extraData.dislocationCaptureRadius(
			_reactants[0], 0);
		rdCl[0][1] = this->_clusterData->extraData.dislocationCaptureRadius(
			_reactants[0], 1);
	}

	if (this->_reactants[1] >= numClusters) {
		auto i = (this->_reactants[1] - numClusters) / 2;
		switch (i) {
		// Void
		case 0:
			r1 = this->_clusterData->bubbleAvRad();
			cl1Reg[Species::V] = {1, 2};
			rdCl[1][0] = 0.79;
			rdCl[1][1] = 1.59;
			break;
		// Perfect V
		case 1:
			r1 = this->_clusterData->perfVAvRad();
			cl1Reg[Species::PerfectV] = {1, 2};
			rdCl[1][0] = 0.79;
			rdCl[1][1] = 1.59;
			break;
		// Faulted V
		case 2:
			r1 = this->_clusterData->faulVAvRad();
			cl1Reg[Species::FaultedV] = {1, 2};
			rdCl[1][0] = 0.79;
			rdCl[1][1] = 1.59;
			break;
		// Perfect I
		case 3:
			r1 = this->_clusterData->perfIAvRad();
			cl1Reg[Species::PerfectI] = {1, 2};
			rdCl[1][0] = 1.85;
			rdCl[1][1] = 1.04;
			break;
		// Faulted I
		case 4:
			r1 = this->_clusterData->faulIAvRad();
			cl1Reg[Species::FaultedI] = {1, 2};
			rdCl[1][0] = 1.85;
			rdCl[1][1] = 1.04;
			break;
		}
	}
	else {
		auto cl1 = this->_clusterData->getCluster(this->_reactants[1]);
		r1 = cl1.getReactionRadius();
		dc1 = cl1.getDiffusionCoefficient(gridIndex);
		cl1Reg = cl1.getRegion();
		rdCl[1][0] = this->_clusterData->extraData.dislocationCaptureRadius(
			_reactants[1], 0);
		rdCl[1][1] = this->_clusterData->extraData.dislocationCaptureRadius(
			_reactants[1], 1);
	}

	return alloy::getRate(cl0Reg, cl1Reg, r0, r1, dc0, dc1, rdCl);
}

KOKKOS_INLINE_FUNCTION
void
AlloyProductionReaction::computeCoefficients()
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
AlloyProductionReaction::computeFlux(
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
		if (comp[Species::V] > 0) {
			// SSBM V case
			if (typeId < 3) {
				// The standard cluster always loses the flux
				Kokkos::atomic_sub(&fluxes[stdClusterId], f);

				// The V size increases (kind of cheating on the ID)
				Kokkos::atomic_add(&fluxes[ssbmId + 1], f * comp[Species::V]);
			}
			// SSBM I case
			else {
				// Special factor to deal with the threshold size
				double gauss = 1.0;

				// Special case where the product is not the single size
				if (this->_products[0] < numClusters) {
					// Only if the large loop has a specific size
					auto avLoop =
						concentrations(ssbmId + 1) / concentrations(ssbmId);
					if (concentrations(ssbmId) == 0.0)
						avLoop = 0.0;

					// Get the product composition
					auto pr =
						this->_clusterData->getCluster(this->_products[0]);
					auto prReg = pr.getRegion();
					Composition prComp(prReg.getOrigin());

					// Target value for the reaction to happen
					double target = prComp[Species::PerfectI] +
						prComp[Species::FaultedI] + comp[Species::V];

					// Gaussian function around it
					double twoSigmaTwo = 0.1;
					gauss = exp(-(avLoop - target) * (avLoop - target) /
								twoSigmaTwo) /
						sqrt(::xolotl::core::pi * twoSigmaTwo);

					// The large bubble concentration decreases
					Kokkos::atomic_sub(&fluxes[ssbmId], f * gauss);
					// The product concentration increases
					Kokkos::atomic_add(&fluxes[this->_products[0]], f * gauss);

					// The I size decreases even more
					Kokkos::atomic_sub(&fluxes[ssbmId + 1],
						f *
							(prComp[Species::PerfectI] +
								prComp[Species::FaultedI]) *
							gauss);
				}

				// In every case

				// The standard cluster always loses the flux
				Kokkos::atomic_sub(&fluxes[stdClusterId], f * gauss);

				// The I size decreases
				Kokkos::atomic_sub(
					&fluxes[ssbmId + 1], f * comp[Species::V] * gauss);
			}
		}

		// Interstitial case
		if (comp[Species::I] > 0) {
			// SSBM V case
			if (typeId < 3) {
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
						prComp[Species::PerfectV] + prComp[Species::FaultedV] +
						comp[Species::I];

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
						f *
							(prComp[Species::V] + prComp[Species::PerfectV] +
								prComp[Species::FaultedV]) *
							gauss);
				}

				// In every case

				// The standard cluster always loses the flux
				Kokkos::atomic_sub(&fluxes[stdClusterId], f * gauss);

				// The V size decreases
				Kokkos::atomic_sub(
					&fluxes[ssbmId + 1], f * comp[Species::I] * gauss);
			}
			// SSBM I case
			else {
				// The standard cluster always loses the flux
				Kokkos::atomic_sub(&fluxes[stdClusterId], f);

				// The I size increases
				Kokkos::atomic_add(&fluxes[ssbmId + 1], f * comp[Species::I]);
			}
		}

		// Helium case
		if (comp[Species::He] > 0) {
			// The standard cluster always loses the flux
			Kokkos::atomic_sub(&fluxes[stdClusterId], f);

			// The He size increases
			Kokkos::atomic_add(&fluxes[ssbmId + 2], f * comp[Species::He]);
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
			auto totalSize = comp1[Species::V] + comp2[Species::V] +
				comp1[Species::PerfectV] + comp2[Species::PerfectV] +
				comp1[Species::FaultedV] + comp2[Species::FaultedV];
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
			auto totalSize = comp1[Species::I] + comp2[Species::I] +
				comp1[Species::PerfectI] + comp2[Species::PerfectI] +
				comp1[Species::FaultedI] + comp2[Species::FaultedI];
			// Both reactants decrease
			Kokkos::atomic_sub(&fluxes[this->_reactants[0]], f);
			Kokkos::atomic_sub(&fluxes[this->_reactants[1]], f);
			// The large bubble increases, as well as average I
			Kokkos::atomic_add(&fluxes[this->_products[0]], f);
			Kokkos::atomic_add(&fluxes[this->_products[0] + 1], f * totalSize);
		}

		// Helium case
		if (orig1.isOnAxis(Species::He) or orig2.isOnAxis(Species::He)) {
			// Compute the total size
			auto totalSize = comp1[Species::He] + comp2[Species::He];
			// Both reactants decrease
			Kokkos::atomic_sub(&fluxes[this->_reactants[0]], f);
			Kokkos::atomic_sub(&fluxes[this->_reactants[1]], f);
			// The large bubble increases, as well as average He
			Kokkos::atomic_add(&fluxes[this->_products[0]], f);
			Kokkos::atomic_add(&fluxes[this->_products[0] + 1], f * totalSize);
		}
	}
}

KOKKOS_INLINE_FUNCTION
void
AlloyProductionReaction::computePartialDerivatives(
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
			if (typeId < 3) {
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
			// SSBM I case
			else {
				// Special factor to deal with the threshold size
				double gauss = 1.0;

				// Special case where the product is not the single size
				if (this->_products[0] < numClusters) {
					// Only if the large void has a specific size
					auto avLoop =
						concentrations(ssbmId + 1) / concentrations(ssbmId);
					if (concentrations(ssbmId) == 0.0)
						avLoop = 0.0;

					// Get the product composition
					auto pr =
						this->_clusterData->getCluster(this->_products[0]);
					auto prReg = pr.getRegion();
					Composition prComp(prReg.getOrigin());

					// Target value for the reaction to happen
					double target = prComp[Species::PerfectI] +
						prComp[Species::FaultedI] + comp[Species::V];

					// Gaussian function around it
					double twoSigmaTwo = 0.1;
					double gauss = exp(-(avLoop - target) * (avLoop - target) /
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

					// The I size decreases even more
					f = this->_coefs(0, 0, 0, 0) * rate *
						(prComp[Species::PerfectI] +
							prComp[Species::FaultedI]) *
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

				// The I size decreases
				f = this->_coefs(0, 0, 0, 0) * rate * comp[Species::V] * gauss;
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
		}
		// Interstitial case
		if (comp[Species::I] > 0) {
			// SSBM V case
			if (typeId < 3) {
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
						prComp[Species::PerfectV] + prComp[Species::FaultedV] +
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
						(prComp[Species::V] + prComp[Species::PerfectV] +
							prComp[Species::FaultedV]) *
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
			// SSBM I case
			else {
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
		}

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
				(comp1[Species::V] + comp2[Species::V] +
					comp1[Species::PerfectV] + comp2[Species::PerfectV] +
					comp1[Species::FaultedV] + comp2[Species::FaultedV]);
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
				(comp1[Species::I] + comp2[Species::I] +
					comp1[Species::PerfectI] + comp2[Species::PerfectI] +
					comp1[Species::FaultedI] + comp2[Species::FaultedI]);
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

	// HeV case
	if (lo[Species::He] > 0) {
		be = prod1.getFormationEnergy() + prod2.getFormationEnergy() -
			cl.getFormationEnergy();
	}
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
	// Standard case
	if (not isLargeBubbleReaction) {
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
	else {
		return this->_clusterData->faulVAv();
	}
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

KOKKOS_INLINE_FUNCTION
void
AlloyTransformReaction::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	// Always call the standard case
	Superclass::computeFlux(concentrations, fluxes, gridIndex);

	if (not isLargeBubbleReaction)
		return;

	auto numClusters = this->_clusterData->numClusters;

	// Switch on V or I type loop
	auto typeId = (_reactant - numClusters) / 2;
	// SSBM V type
	if (typeId < 3) {
		// The size of the Perfect loop increases by the size of the Faulted
		// loop
		auto fSize = this->_clusterData->faulVAv();
		Kokkos::atomic_add(&fluxes(this->_clusterData->perfVAvId()),
			this->_rate(gridIndex) * concentrations(_reactant) * fSize);
	}
	// SSBM I type
	else {
		// The size of the Perfect loop increases by the size of the Faulted
		// loop
		auto fSize = this->_clusterData->faulIAv();
		Kokkos::atomic_add(&fluxes(this->_clusterData->perfIAvId()),
			this->_rate(gridIndex) * concentrations(_reactant) * fSize);
	}
}

KOKKOS_INLINE_FUNCTION
void
AlloyTransformReaction::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	// Always call the standard case
	Superclass::computePartialDerivatives(concentrations, values, gridIndex);

	if (not isLargeBubbleReaction)
		return;

	auto numClusters = this->_clusterData->numClusters;

	// Switch on V or I type loop
	auto typeId = (_reactant - numClusters) / 2;

	double fSize = 0.0;
	// SSBM V type
	if (typeId < 3) {
		fSize = this->_clusterData->faulVAv();
	}
	// SSBM I type
	else {
		fSize = this->_clusterData->faulIAv();
	}
	// The size of the Perfect loop increases by the size of the Faulted loop
	Kokkos::atomic_add(
		&values(_connEntries[1][1][0][0]), this->_rate(gridIndex) * fSize);
}
} // namespace network
} // namespace core
} // namespace xolotl
