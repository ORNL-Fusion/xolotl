#pragma once

#include <xolotl/core/network/detail/impl/SinkReactionGenerator.tpp>
#include <xolotl/core/network/detail/impl/TransformReactionGenerator.tpp>
#include <xolotl/core/network/impl/FeCrClusterGenerator.tpp>
#include <xolotl/core/network/impl/FeCrReaction.tpp>
#include <xolotl/core/network/impl/ReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
void
FeCrReactionNetwork::initializeExtraClusterData(
	const options::IOptions& options)
{
	this->_clusterData.h_view().extraData.initialize(
		this->_clusterData.h_view().numClusters);
	this->copyClusterDataView();
	this->invalidateDataMirror();
}

void
FeCrReactionNetwork::computeFluxesPreProcess(ConcentrationsView concentrations,
	FluxesView fluxes, IndexType gridIndex, double surfaceDepth, double spacing)
{
	auto data = this->_clusterData.h_view();
	auto diffusionFactor = this->getClusterDataMirror().diffusionFactor;
	auto sigma = create_mirror_view(data.extraData.netSigma);
	for (auto i = 0; i < this->_numClusters; i++) {
		if (diffusionFactor(i) == 0.0) {
			sigma(i) = 0.0;
			continue;
		}
		sigma(i) = this->getNetSigma(concentrations, i, gridIndex);
	}
	deep_copy(data.extraData.netSigma, sigma);
	this->updateReactionRates(0.0);
}

void
FeCrReactionNetwork::computePartialsPreProcess(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex, double surfaceDepth, double spacing)
{
	auto data = this->_clusterData.h_view();
	auto diffusionFactor = this->getClusterDataMirror().diffusionFactor;
	auto sigma = create_mirror_view(data.extraData.netSigma);
	for (auto i = 0; i < this->_numClusters; i++) {
		if (diffusionFactor(i) == 0.0) {
			sigma(i) = 0.0;
			continue;
		}
		sigma(i) = this->getNetSigma(concentrations, i, gridIndex);
	}
	deep_copy(data.extraData.netSigma, sigma);
	this->updateReactionRates(0.0);
}

namespace detail
{
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
FeCrReactionGenerator::operator()(IndexType i, IndexType j, TTag tag) const
{
	// Transform happens between sessile clusters
	addTransforms(i, j, tag);

	// Check the diffusion factors
	auto diffusionFactor = this->_clusterData.diffusionFactor;
	if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
		return;
	}

	using Species = typename Network::Species;
	using Composition = typename Network::Composition;
	using AmountType = typename Network::AmountType;

	if (i == j) {
		addSinks(i, tag);
	}

	auto& subpaving = this->getSubpaving();
	auto previousIndex = subpaving.invalidIndex();
	auto numClusters = this->getNumberOfClusters();

	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	// void + void → void
	if (lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::V)) {
		// Find the grouped one
		auto groupedReg = lo1[Species::V] > lo2[Species::V] ? cl1Reg : cl2Reg;
		double singleSize = lo1[Species::V] > lo2[Species::V] ?
			lo2[Species::V] :
			lo1[Species::V];

		for (auto k : makeIntervalRange(groupedReg[Species::V])) {
			// Compute the composition of the new cluster
			auto size = singleSize + k;
			// Find the corresponding cluster
			Composition comp = Composition::zero();
			comp[Species::V] = size;
			auto vProdId = subpaving.findTileId(comp);
			if (vProdId != subpaving.invalidIndex() &&
				vProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, vProdId});
				previousIndex = vProdId;
				if (lo1[Species::V] == 1 || lo2[Species::V] == 1) {
					this->addDissociationReaction(tag, {vProdId, i, j});
				}
			}
		}
		return;
	}

	// void + int → void | int | recombine
	if ((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::I))) {
		// Find out which one is which
		auto vReg = lo1.isOnAxis(Species::V) ? cl1Reg : cl2Reg;
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];

		for (auto k : makeIntervalRange(vReg[Species::V])) {
			// Compute the product size
			int prodSize = k - iSize;
			// 3 cases
			if (prodSize > 0) {
				// Looking for V cluster
				Composition comp = Composition::zero();
				comp[Species::V] = prodSize;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() &&
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					previousIndex = vProdId;
					// No dissociation
				}
			}
			else if (prodSize < 0) {
				// Looking for I cluster
				Composition comp = Composition::zero();
				comp[Species::I] = -prodSize;
				auto iProdId = subpaving.findTileId(comp);
				if (iProdId != subpaving.invalidIndex()) {
					this->addProductionReaction(tag, {i, j, iProdId});
					// No dissociation
				}
			}
			else {
				// No product
				this->addProductionReaction(tag, {i, j});
			}
		}
		return;
	}

	// void + trapped → void | trapped | complex | trap
	if ((lo1[Species::Trapped] > 0 && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2[Species::Trapped] > 0)) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		auto tSize = lo1.isOnAxis(Species::V) ? lo2[Species::Trapped] :
												lo1[Species::Trapped];
		// Compute the product size
		int prodSize = tSize - vSize;
		if (prodSize > 0) {
			//  Find the corresponding cluster
			Composition comp = Composition::zero();
			comp[Species::Trapped] = prodSize;
			comp[Species::Trap] = 1;
			auto tProdId = subpaving.findTileId(comp);
			if (tProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, tProdId});
				// No dissociation
			}
			comp[Species::Trapped] = 0;
			comp[Species::Complex] = prodSize;
			auto cProdId = subpaving.findTileId(comp);
			if (cProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, cProdId});
				// No dissociation
			}
		}
		else if (prodSize < 0) {
			// Looking for V cluster
			Composition comp = Composition::zero();
			comp[Species::V] = -prodSize;
			auto vProdId = subpaving.findTileId(comp);
			if (vProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, vProdId});
				// No dissociation
			}
		}
		else {
			// Trap is the product
			Composition comp = Composition::zero();
			comp[Species::Trap] = 1;
			auto tProdId = subpaving.findTileId(comp);
			if (tProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, tProdId});
				// No dissociation
			}
		}
		return;
	}

	// void + free → void | free | int | recombine
	if ((lo1.isOnAxis(Species::Free) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::Free))) {
		// Find out which one is which
		auto vReg = lo1.isOnAxis(Species::V) ? cl1Reg : cl2Reg;
		auto fReg = lo1.isOnAxis(Species::V) ? cl2Reg : cl1Reg;

		// Compute the minimum and maximum size (on the I axis)
		int minSize = fReg[Species::Free].begin() - vReg[Species::V].end() + 1;
		int maxSize = fReg[Species::Free].end() - 1 - vReg[Species::V].begin();

		// Recombine
		if (minSize <= 0 and maxSize >= 0) {
			this->addProductionReaction(tag, {i, j});
		}

		// Look for appropriate products
		for (IndexType k = 0; k < numClusters; ++k) {
			// Get the composition
			const auto& prodReg = this->getCluster(k).getRegion();
			Composition loProd = prodReg.getOrigin();
			bool isGood = true;

			// Void
			if (minSize < 0) {
				if (!loProd.isOnAxis(Species::V))
					continue;
				if (prodReg[Species::V].begin() > -minSize) {
					isGood = false;
					continue;
				}
				if (prodReg[Species::V].end() - 1 < -util::min(0, maxSize)) {
					isGood = false;
					continue;
				}

				if (isGood) {
					this->addProductionReaction(tag, {i, j, k});
				}
			}
			// Free and Int
			if (maxSize > 0) {
				if (!loProd.isOnAxis(Species::I) and
					!loProd.isOnAxis(Species::Free))
					continue;
				if (prodReg[Species::I].begin() +
						prodReg[Species::Free].begin() >
					maxSize) {
					isGood = false;
					continue;
				}
				if (prodReg[Species::I].end() + prodReg[Species::Free].end() -
						2 <
					util::max(0, minSize)) {
					isGood = false;
					continue;
				}

				if (isGood) {
					this->addProductionReaction(tag, {i, j, k});
				}
			}
		}
		return;
	}

	// void + loop → loop | complex | trap
	if ((lo1[Species::Loop] > 0 && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2[Species::Loop] > 0)) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		auto lReg = lo1.isOnAxis(Species::V) ? cl2Reg : cl1Reg;

		for (auto k : makeIntervalRange(lReg[Species::Loop])) {
			// Compute the product size
			int prodSize = k - vSize;
			if (prodSize > 0) {
				//  Find the corresponding cluster
				Composition comp = Composition::zero();
				comp[Species::Loop] = prodSize;
				comp[Species::Trap] = 1;
				auto lProdId = subpaving.findTileId(comp);
				if (lProdId != subpaving.invalidIndex() &&
					lProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, lProdId});
					previousIndex = lProdId;
					// No dissociation
				}
				comp[Species::Loop] = 0;
				comp[Species::Complex] = prodSize;
				auto cProdId = subpaving.findTileId(comp);
				if (cProdId != subpaving.invalidIndex() &&
					cProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, cProdId});
					previousIndex = cProdId;
					// No dissociation
				}
			}
			else if (prodSize == 0) {
				// Trap is the product
				Composition comp = Composition::zero();
				comp[Species::Trap] = 1;
				auto tProdId = subpaving.findTileId(comp);
				if (tProdId != subpaving.invalidIndex() &&
					tProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, tProdId});
					previousIndex = tProdId;
					// No dissociation
				}
			}
		}
		return;
	}

	// void + junction → junction | trapped | complex | trap
	if ((lo1[Species::Junction] > 0 && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2[Species::Junction] > 0)) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		auto jReg = lo1.isOnAxis(Species::V) ? cl2Reg : cl1Reg;

		for (auto k : makeIntervalRange(jReg[Species::Junction])) {
			// Compute the product size
			int prodSize = k - vSize;
			if (prodSize > 0) {
				//  Find the corresponding cluster
				Composition comp = Composition::zero();
				comp[Species::Junction] = prodSize;
				comp[Species::Trap] = 1;
				auto jProdId = subpaving.findTileId(comp);
				if (jProdId != subpaving.invalidIndex() &&
					jProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, jProdId});
					previousIndex = jProdId;
					// No dissociation
					return;
				}
				comp[Species::Junction] = 0;
				comp[Species::Trapped] = prodSize;
				auto tProdId = subpaving.findTileId(comp);
				if (tProdId != subpaving.invalidIndex() &&
					tProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, tProdId});
					previousIndex = tProdId;
					// No dissociation
					return;
				}
				comp[Species::Trapped] = 0;
				comp[Species::Complex] = prodSize;
				auto cProdId = subpaving.findTileId(comp);
				if (cProdId != subpaving.invalidIndex() &&
					cProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, cProdId});
					previousIndex = cProdId;
					// No dissociation
				}
			}
			// This never actually happens because only small V interact and
			// Junctions are too big
			else if (prodSize == 0) {
				// Trap is the product
				Composition comp = Composition::zero();
				comp[Species::Trap] = 1;
				auto tProdId = subpaving.findTileId(comp);
				if (tProdId != subpaving.invalidIndex() &&
					tProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, tProdId});
					previousIndex = tProdId;
					// No dissociation
				}
			}
		}
		return;
	}

	// int + int → int | free
	if (lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::I)) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::I] + lo2[Species::I];
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::I] = size;
		auto iProdId = subpaving.findTileId(comp);
		if (iProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, iProdId});
			if (lo1[Species::I] == 1 || lo2[Species::I] == 1) {
				this->addDissociationReaction(tag, {iProdId, i, j});
			}
		}
		comp[Species::I] = 0;
		comp[Species::Free] = size;
		auto fProdId = subpaving.findTileId(comp);
		if (fProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, fProdId});
			if (lo1[Species::I] == 1 || lo2[Species::I] == 1) {
				this->addDissociationReaction(tag, {fProdId, i, j});
			}
		}
		return;
	}

	// int + trapped → trapped
	if ((lo1[Species::Trapped] > 0 && lo2.isOnAxis(Species::I)) ||
		(lo1.isOnAxis(Species::I) && lo2[Species::Trapped] > 0)) {
		// Find out which one is which
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];
		auto tReg = lo1.isOnAxis(Species::I) ? cl2Reg : cl1Reg;

		for (auto k : makeIntervalRange(tReg[Species::Trapped])) {
			// Compute the composition of the new cluster
			auto size = iSize + k;
			// Find the corresponding cluster
			Composition comp = Composition::zero();
			comp[Species::Trapped] = size;
			comp[Species::Trap] = 1;
			auto tProdId = subpaving.findTileId(comp);
			if (tProdId != subpaving.invalidIndex() &&
				tProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, tProdId});
				previousIndex = tProdId;
				if (iSize == 1) {
					this->addDissociationReaction(tag, {tProdId, i, j});
				}
			}
		}
		return;
	}

	// int + free → free
	if ((lo1.isOnAxis(Species::Free) && lo2.isOnAxis(Species::I)) ||
		(lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::Free))) {
		// Find out which one is which
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];
		auto fReg = lo1.isOnAxis(Species::I) ? cl2Reg : cl1Reg;

		for (auto k : makeIntervalRange(fReg[Species::Free])) {
			// Compute the composition of the new cluster
			auto size = iSize + k;
			// Find the corresponding cluster
			Composition comp = Composition::zero();
			comp[Species::Free] = size;
			auto fProdId = subpaving.findTileId(comp);
			if (fProdId != subpaving.invalidIndex() &&
				fProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, fProdId});
				previousIndex = fProdId;
				if (lo1[Species::I] == 1 || lo2[Species::I] == 1) {
					this->addDissociationReaction(tag, {fProdId, i, j});
				}
			}
		}
		return;
	}

	// int + loop → loop
	if ((lo1[Species::Loop] > 0 && lo2.isOnAxis(Species::I)) ||
		(lo1.isOnAxis(Species::I) && lo2[Species::Loop] > 0)) {
		// Find out which one is which
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];
		auto lReg = lo1.isOnAxis(Species::I) ? cl2Reg : cl1Reg;

		for (auto k : makeIntervalRange(lReg[Species::Loop])) {
			// Compute the composition of the new cluster
			auto size = iSize + k;
			// Find the corresponding cluster
			Composition comp = Composition::zero();
			comp[Species::Loop] = size;
			comp[Species::Trap] = 1;
			auto lProdId = subpaving.findTileId(comp);
			if (lProdId != subpaving.invalidIndex() &&
				lProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, lProdId});
				previousIndex = lProdId;
				// No dissociation
			}
		}
		return;
	}

	// int + junction → junction
	if ((lo1[Species::Junction] > 0 && lo2.isOnAxis(Species::I)) ||
		(lo1.isOnAxis(Species::I) && lo2[Species::Junction] > 0)) {
		// Find out which one is which
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];
		auto jReg = lo1.isOnAxis(Species::I) ? cl2Reg : cl1Reg;

		for (auto k : makeIntervalRange(jReg[Species::Junction])) {
			// Compute the composition of the new cluster
			auto size = iSize + k;
			// Find the corresponding cluster
			Composition comp = Composition::zero();
			comp[Species::Junction] = size;
			comp[Species::Trap] = 1;
			auto jProdId = subpaving.findTileId(comp);
			if (jProdId != subpaving.invalidIndex() &&
				jProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, jProdId});
				previousIndex = jProdId;
				// No dissociation
			}
		}
		return;
	}

	// free + trapped → junction | trapped
	if ((lo1[Species::Trapped] > 0 && lo2.isOnAxis(Species::Free)) ||
		(lo1.isOnAxis(Species::Free) && lo2[Species::Trapped] > 0)) {
		auto previousJunctionId = subpaving.invalidIndex();
		auto previousTrappedId = subpaving.invalidIndex();

		// Find out which one is which
		auto fReg = lo1.isOnAxis(Species::Free) ? cl1Reg : cl2Reg;
		auto tReg = lo1.isOnAxis(Species::Free) ? cl2Reg : cl1Reg;

		// Compute the minimum and maximum size (on the I axis)
		int minSize =
			fReg[Species::Free].begin() + tReg[Species::Trapped].begin();
		int maxSize =
			fReg[Species::Free].end() + tReg[Species::Trapped].end() - 2;

		double ratio = 2.0 *
			fabs((double)fReg[Species::Free].begin() -
				(double)tReg[Species::Trapped].begin()) /
			((double)fReg[Species::Free].begin() +
				(double)tReg[Species::Trapped].begin());

		for (auto size = minSize; size <= maxSize; size++) {
			// Junction case
			if (ratio < 0.5) {
				Composition comp = Composition::zero();
				comp[Species::Junction] = size;
				comp[Species::Trap] = 1;
				auto jProdId = subpaving.findTileId(comp);
				if (jProdId != subpaving.invalidIndex() and
					jProdId != previousJunctionId) {
					this->addProductionReaction(tag, {i, j, jProdId});
					previousJunctionId = jProdId;
					// No dissociation
				}
			}
			// Trapped case
			if (ratio >= 0.5) {
				Composition comp = Composition::zero();
				comp[Species::Trapped] = size;
				comp[Species::Trap] = 1;
				auto tProdId = subpaving.findTileId(comp);
				if (tProdId != subpaving.invalidIndex() and
					tProdId != previousTrappedId) {
					this->addProductionReaction(tag, {i, j, tProdId});
					previousTrappedId = tProdId;
					// No dissociation
				}
			}
		}
		return;
	}

	// free + trap → trapped
	if ((lo1.isOnAxis(Species::Free) && lo2.isOnAxis(Species::Trap)) ||
		(lo1.isOnAxis(Species::Trap) && lo2.isOnAxis(Species::Free))) {
		// Find out which one is which
		auto fReg = lo1.isOnAxis(Species::Free) ? cl1Reg : cl2Reg;

		for (auto k : makeIntervalRange(fReg[Species::Free])) {
			// Compute the composition of the new cluster
			auto size = k;
			// Find the corresponding cluster
			Composition comp = Composition::zero();
			comp[Species::Trapped] = size;
			comp[Species::Trap] = 1;
			auto tProdId = subpaving.findTileId(comp);
			if (tProdId != subpaving.invalidIndex() &&
				tProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, tProdId});
				previousIndex = tProdId;
				this->addDissociationReaction(tag, {tProdId, i, j});

				// Special dissociation
				comp[Species::Trapped] = 0;
				comp[Species::Junction] = size;
				auto jProdId = subpaving.findTileId(comp);
				if (jProdId != subpaving.invalidIndex()) {
					this->addDissociationReaction(tag, {jProdId, i, j});
				}
			}
		}
		return;
	}

	// free + loop → junction | trapped | loop
	if ((lo1[Species::Loop] > 0 && lo2.isOnAxis(Species::Free)) ||
		(lo1.isOnAxis(Species::Free) && lo2[Species::Loop] > 0)) {
		auto previousJunctionId = subpaving.invalidIndex();
		auto previousLoopId = subpaving.invalidIndex();
		auto previousTrappedId = subpaving.invalidIndex();

		// Find out which one is which
		auto fReg = lo1.isOnAxis(Species::Free) ? cl1Reg : cl2Reg;
		auto lReg = lo1.isOnAxis(Species::Free) ? cl2Reg : cl1Reg;

		// Compute the minimum and maximum size (on the I axis)
		int minSize = fReg[Species::Free].begin() + lReg[Species::Loop].begin();
		int maxSize = fReg[Species::Free].end() + lReg[Species::Loop].end() - 2;

		double ratio = 2.0 *
			fabs((double)fReg[Species::Free].begin() -
				(double)lReg[Species::Loop].begin()) /
			((double)fReg[Species::Free].begin() +
				(double)lReg[Species::Loop].begin());

		for (auto size = minSize; size <= maxSize; size++) {
			// Junction case
			if (ratio < 0.5) {
				Composition comp = Composition::zero();
				comp[Species::Junction] = size;
				comp[Species::Trap] = 1;
				auto jProdId = subpaving.findTileId(comp);
				if (jProdId != subpaving.invalidIndex() and
					jProdId != previousJunctionId) {
					this->addProductionReaction(tag, {i, j, jProdId});
					previousJunctionId = jProdId;
					// No dissociation
				}
			}
			// Trapped or Loop case
			if (ratio > 0.5) {
				// Loop
				if (fReg[Species::Free].begin() < lReg[Species::Loop].begin()) {
					Composition comp = Composition::zero();
					comp[Species::Loop] = size;
					comp[Species::Trap] = 1;
					auto lProdId = subpaving.findTileId(comp);
					if (lProdId != subpaving.invalidIndex() and
						lProdId != previousLoopId) {
						this->addProductionReaction(tag, {i, j, lProdId});
						previousLoopId = lProdId;
						// No dissociation
					}
				}
				// Trapped
				else {
					Composition comp = Composition::zero();
					comp[Species::Trapped] = size;
					comp[Species::Trap] = 1;
					auto tProdId = subpaving.findTileId(comp);
					if (tProdId != subpaving.invalidIndex() and
						tProdId != previousTrappedId) {
						this->addProductionReaction(tag, {i, j, tProdId});
						previousTrappedId = tProdId;
						// No dissociation
					}
				}
			}
		}
		return;
	}

	// free + junction → junction
	if ((lo1[Species::Junction] > 0 && lo2.isOnAxis(Species::Free)) ||
		(lo1.isOnAxis(Species::Free) && lo2[Species::Junction] > 0)) {
		// Find out which one is which
		auto fReg = lo1.isOnAxis(Species::Free) ? cl1Reg : cl2Reg;
		auto jReg = lo1.isOnAxis(Species::Free) ? cl2Reg : cl1Reg;

		// Compute the minimum and maximum size (on the I axis)
		int minSize =
			fReg[Species::Free].begin() + jReg[Species::Junction].begin();
		int maxSize =
			fReg[Species::Free].end() + jReg[Species::Junction].end() - 2;

		// Look for appropriate products
		for (IndexType k = 0; k < numClusters; ++k) {
			// Get the composition
			const auto& prodReg = this->getCluster(k).getRegion();
			Composition loProd = prodReg.getOrigin();
			bool isGood = true;

			// Junction
			if (loProd[Species::Junction] == 0)
				continue;
			if (prodReg[Species::Junction].begin() > maxSize) {
				isGood = false;
				continue;
			}
			if (prodReg[Species::Junction].end() - 1 < minSize) {
				isGood = false;
				continue;
			}

			if (isGood) {
				this->addProductionReaction(tag, {i, j, k});
			}
		}
		return;
	}

	// Disso: complex → trap + int
	if ((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::Trap)) ||
		(lo1.isOnAxis(Species::Trap) && lo2.isOnAxis(Species::I))) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::I] +
			lo2[Species::I]; // The other axis should be 0 so it should
							 // work to add up everything
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::Complex] = size;
		comp[Species::Trap] = 1;
		auto cProdId = subpaving.findTileId(comp);
		if (cProdId != subpaving.invalidIndex()) {
			this->addDissociationReaction(tag, {cProdId, i, j});
		}
		return;
	}

	return;
}

template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
FeCrReactionGenerator::addSinks(IndexType i, TTag tag) const
{
	using Species = typename Network::Species;
	using Composition = typename Network::Composition;

	const auto& clReg = this->getCluster(i).getRegion();
	Composition lo = clReg.getOrigin();

	// I
	if (clReg.isSimplex() && lo.isOnAxis(Species::I)) {
		this->addSinkReaction(tag, {i, Network::invalidIndex()});
	}

	// Void
	if (clReg.isSimplex() && lo.isOnAxis(Species::V) && lo[Species::V] < 5) {
		this->addSinkReaction(tag, {i, Network::invalidIndex()});
	}

	// Free
	if (lo.isOnAxis(Species::Free)) {
		this->addSinkReaction(tag, {i, Network::invalidIndex()});
	}
}

template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
FeCrReactionGenerator::addTransforms(IndexType i, IndexType j, TTag tag) const
{
	using Species = typename Network::Species;
	using Composition = typename Network::Composition;

	const auto& cl1Reg = this->getCluster(i).getRegion();
	Composition lo1 = cl1Reg.getOrigin();

	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo2 = cl2Reg.getOrigin();

	// The first reactant need to be Junction
	if (lo1[Species::Junction] == 0 and lo2[Species::Junction] == 0)
		return;

	// Which is Junction ?
	auto junction = lo1[Species::Junction] > 0 ? i : j;
	auto junctionReg = lo1[Species::Junction] > 0 ? lo1 : lo2;
	auto other = lo1[Species::Junction] > 0 ? j : i;
	auto otherReg = lo1[Species::Junction] > 0 ? lo2 : lo1;

	// The second on should be Trapped or Loop
	if (otherReg[Species::Trapped] == 0 && otherReg[Species::Loop] == 0)
		return;

	// They need to be the same size
	// TODO: this only works is the grouping is the same in every direction
	if (otherReg[Species::Trapped] == junctionReg[Species::Junction] ||
		otherReg[Species::Loop] == junctionReg[Species::Junction]) {
		this->addTransformReaction(tag, {junction, other});
	}
}

inline ReactionCollection<FeCrReactionGenerator::Network>
FeCrReactionGenerator::getReactionCollection() const
{
	ReactionCollection<Network> ret(this->_clusterData.gridSize,
		this->getProductionReactions(), this->getDissociationReactions(),
		this->getSinkReactions(), this->getTransformReactions());
	return ret;
}
} // namespace detail

inline detail::FeCrReactionGenerator
FeCrReactionNetwork::getReactionGenerator() const noexcept
{
	return detail::FeCrReactionGenerator{*this};
}
} // namespace network
} // namespace core
} // namespace xolotl
