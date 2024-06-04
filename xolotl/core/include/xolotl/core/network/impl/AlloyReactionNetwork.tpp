#pragma once

#include <xolotl/core/network/detail/impl/SinkReactionGenerator.tpp>
#include <xolotl/core/network/detail/impl/TransformReactionGenerator.tpp>
#include <xolotl/core/network/impl/AlloyClusterGenerator.tpp>
#include <xolotl/core/network/impl/AlloyReaction.tpp>
#include <xolotl/core/network/impl/ReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
AlloyReactionGenerator::operator()(IndexType i, IndexType j, TTag tag) const
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

	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	// vac + vac = vac | faultedV
	if (lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::V) &&
		diffusionFactor(i) > 0.0 && diffusionFactor(j) > 0.0) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::V] + lo2[Species::V];
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::V] = size;
		auto vProdId = subpaving.findTileId(comp);
		if (vProdId != subpaving.invalidIndex()) {
			// Check the diffusivity to distinguish V from Voids
			if (diffusionFactor(vProdId) > 0.0) {
				this->addProductionReaction(tag, {i, j, vProdId});
				if (lo1[Species::V] == 1 || lo2[Species::V] == 1) {
					this->addDissociationReaction(tag, {vProdId, i, j});
				}
			}
		}
		comp[Species::V] = 0;
		comp[Species::FaultedV] = size;
		auto fProdId = subpaving.findTileId(comp);
		if (fProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, fProdId});
			if (lo1[Species::V] == 1 || lo2[Species::V] == 1) {
				this->addDissociationReaction(tag, {fProdId, i, j});
			}
		}
		return;
	}

	// vac + void = void
	if (lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::V) &&
		(diffusionFactor(i) > 0.0 || diffusionFactor(j) > 0.0)) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::V] + lo2[Species::V];
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::V] = size;
		auto vProdId = subpaving.findTileId(comp);
		if (vProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, vProdId});
			if (lo1[Species::V] == 1 || lo2[Species::V] == 1) {
				this->addDissociationReaction(tag, {vProdId, i, j});
			}
		}
		return;
	}

	// vac + perfectV = perfectV
	if ((lo1.isOnAxis(Species::PerfectV) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::PerfectV))) {
		// PerfectV can be grouped
		auto minSize = lo1[Species::V] + lo2[Species::V] +
			lo1[Species::PerfectV] + lo2[Species::PerfectV];
		auto maxSize = lo1[Species::V] + lo2[Species::V] +
			hi1[Species::PerfectV] + hi2[Species::PerfectV] -
			2; // The other axis should be 0 so it should
			   // work to add up everything
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::PerfectV] = k;
			auto pProdId = subpaving.findTileId(comp);
			if (pProdId != subpaving.invalidIndex() &&
				pProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, pProdId});
				if (lo1[Species::V] == 1 || lo2[Species::V] == 1) {
					this->addDissociationReaction(tag, {pProdId, i, j});
				}
				previousIndex = pProdId;
			}
		}
		return;
	}

	// vac + faultedV = faultedV
	if ((lo1.isOnAxis(Species::FaultedV) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::FaultedV))) {
		// Faulted can be grouped
		auto minSize = lo1[Species::V] + lo2[Species::V] +
			lo1[Species::FaultedV] + lo2[Species::FaultedV];
		auto maxSize = lo1[Species::V] + lo2[Species::V] +
			hi1[Species::FaultedV] + hi2[Species::FaultedV] -
			2; // The other axis should be 0 so it should work
			   // to add up everything
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::FaultedV] = k;
			auto fProdId = subpaving.findTileId(comp);
			if (fProdId != subpaving.invalidIndex() &&
				fProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, fProdId});
				if (lo1[Species::V] == 1 || lo2[Species::V] == 1) {
					this->addDissociationReaction(tag, {fProdId, i, j});
				}
				previousIndex = fProdId;
			}
		}
		return;
	}

	// vac + int = vac | int | recombine
	if (((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
			(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::I)))) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];
		// Compute the product size
		int prodSize = vSize - iSize;
		// 3 cases
		if (prodSize > 0) {
			// Looking for V cluster
			Composition comp = Composition::zero();
			comp[Species::V] = prodSize;
			auto vProdId = subpaving.findTileId(comp);
			if (vProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, vProdId});
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
		return;
	}

	// vac + perfectI = perfectI | int | recombine | vac
	if (((lo1.isOnAxis(Species::PerfectI) && lo2.isOnAxis(Species::V)) ||
			(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::PerfectI)))) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		// Perfect can be grouped
		auto pReg = lo1.isOnAxis(Species::PerfectI) ? cl1Reg : cl2Reg;
		for (auto k : makeIntervalRange(pReg[Species::PerfectI])) {
			// Compute the product size
			int prodSize = k - vSize;
			if (prodSize > 0) {
				//  Find the corresponding cluster
				Composition comp = Composition::zero();
				comp[Species::I] = prodSize;
				auto iProdId = subpaving.findTileId(comp);
				if (iProdId != subpaving.invalidIndex() &&
					iProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, iProdId});
					// No dissociation
					previousIndex = iProdId;
				}
				comp[Species::I] = 0;
				comp[Species::PerfectI] = prodSize;
				auto pProdId = subpaving.findTileId(comp);
				if (pProdId != subpaving.invalidIndex() &&
					pProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, pProdId});
					// No dissociation
					previousIndex = pProdId;
				}
			}
			else if (prodSize < 0) {
				// Looking for V cluster
				Composition comp = Composition::zero();
				comp[Species::V] = -prodSize;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() &&
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					// No dissociation
					previousIndex = vProdId;
				}
			}
			else {
				// No product
				this->addProductionReaction(tag, {i, j});
			}
		}
		return;
	}

	// vac + faultedI = faultedI | int | recombine | vac
	if (((lo1.isOnAxis(Species::FaultedI) && lo2.isOnAxis(Species::V)) ||
			(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::FaultedI)))) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		// Frank can be grouped
		auto fReg = lo1.isOnAxis(Species::FaultedI) ? cl1Reg : cl2Reg;
		for (auto k : makeIntervalRange(fReg[Species::FaultedI])) {
			// Compute the product size
			int prodSize = k - vSize;
			if (prodSize > 0) {
				//  Find the corresponding cluster
				Composition comp = Composition::zero();
				comp[Species::I] = prodSize;
				auto iProdId = subpaving.findTileId(comp);
				if (iProdId != subpaving.invalidIndex() &&
					iProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, iProdId});
					previousIndex = iProdId;
					// No dissociation
				}
				comp[Species::I] = 0;
				comp[Species::FaultedI] = prodSize;
				auto fProdId = subpaving.findTileId(comp);
				if (fProdId != subpaving.invalidIndex() &&
					fProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, fProdId});
					previousIndex = fProdId;
					// No dissociation
				}
			}
			else if (prodSize < 0) {
				// Looking for V cluster
				Composition comp = Composition::zero();
				comp[Species::V] = -prodSize;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() &&
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					previousIndex = vProdId;
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

	// int + int = int | faultedI and perfectI
	if (lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::I)) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::I] + lo2[Species::I];
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::I] = size;
		auto iProdId = subpaving.findTileId(comp);
		if (iProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, iProdId});
			//			if (lo1[Species::I] == 1 || lo2[Species::I] == 1) {
			//				this->addDissociationReaction(tag, {iProdId, i, j});
			//			}
		}
		comp[Species::I] = 0;
		comp[Species::FaultedI] = size;
		auto fProdId = subpaving.findTileId(comp);
		if (fProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, fProdId});
			// No dissociation
		}
		comp[Species::FaultedI] = 0;
		comp[Species::PerfectI] = size;
		auto pProdId = subpaving.findTileId(comp);
		if (pProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, pProdId});
			// No dissociation
		}
		return;
	}

	// int + faultedI = faultedI
	if ((lo1.isOnAxis(Species::FaultedI) && lo2.isOnAxis(Species::I)) ||
		(lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::FaultedI))) {
		// Frank can be grouped
		auto minSize = lo1[Species::I] + lo2[Species::I] +
			lo1[Species::FaultedI] + lo2[Species::FaultedI];
		auto maxSize = lo1[Species::I] + lo2[Species::I] +
			hi1[Species::FaultedI] + hi2[Species::FaultedI] -
			2; // The other axis should be 0 so it should
			   // work to add up everything
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::FaultedI] = k;
			auto fProdId = subpaving.findTileId(comp);
			if (fProdId != subpaving.invalidIndex() &&
				fProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, fProdId});
				previousIndex = fProdId;
				// No dissociation
			}
		}
		return;
	}

	// int + perfectI = perfectI
	if ((lo1.isOnAxis(Species::PerfectI) && lo2.isOnAxis(Species::I)) ||
		(lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::PerfectI))) {
		// Perfect can be grouped
		auto minSize = lo1[Species::I] + lo2[Species::I] +
			lo1[Species::PerfectI] + lo2[Species::PerfectI];
		auto maxSize = lo1[Species::I] + lo2[Species::I] +
			hi1[Species::PerfectI] + hi2[Species::PerfectI] -
			2; // The other axis should be 0 so it should
			   // work to add up everything
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::PerfectI] = k;
			auto pProdId = subpaving.findTileId(comp);
			if (pProdId != subpaving.invalidIndex() &&
				pProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, pProdId});
				previousIndex = pProdId;
				// No dissociation
			}
		}
		return;
	}

	// int + faultedV = faultedV | vac
	if (((lo1.isOnAxis(Species::FaultedV) && lo2.isOnAxis(Species::I)) ||
			(lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::FaultedV)))) {
		// Find out which one is which
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];
		// Faulted can be grouped
		auto fReg = lo1.isOnAxis(Species::FaultedV) ? cl1Reg : cl2Reg;
		for (auto k : makeIntervalRange(fReg[Species::FaultedV])) {
			// Compute the product size
			int prodSize = k - iSize;
			if (prodSize > 0) {
				//  Find the corresponding cluster
				Composition comp = Composition::zero();
				comp[Species::V] = prodSize;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() &&
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					previousIndex = vProdId;
					// No dissociation
				}
				comp[Species::V] = 0;
				comp[Species::FaultedV] = prodSize;
				auto fProdId = subpaving.findTileId(comp);
				if (fProdId != subpaving.invalidIndex() &&
					fProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, fProdId});
					previousIndex = fProdId;
					// No dissociation
				}
			}
		}
		return;
	}

	// int + perfectV = perfectV | vac
	if (((lo1.isOnAxis(Species::PerfectV) && lo2.isOnAxis(Species::I)) ||
			(lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::PerfectV)))) {
		// Find out which one is which
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];
		// PerfectV can be grouped
		auto vReg = lo1.isOnAxis(Species::PerfectV) ? cl1Reg : cl2Reg;
		for (auto k : makeIntervalRange(vReg[Species::PerfectV])) {
			// Compute the product size
			int prodSize = k - iSize;
			if (prodSize > 0) {
				//  Find the corresponding cluster
				Composition comp = Composition::zero();
				comp[Species::V] = prodSize;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() &&
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					previousIndex = vProdId;
					// No dissociation
				}
				comp[Species::V] = 0;
				comp[Species::PerfectV] = prodSize;
				auto pProdId = subpaving.findTileId(comp);
				if (pProdId != subpaving.invalidIndex() &&
					pProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, pProdId});
					//					if (lo1[Species::I] == 1 ||
					// lo2[Species::I]
					//== 1) {
					// this->addDissociationReaction(tag, {pProdId, i, j});
					//					}
					previousIndex = pProdId;
				}
			}
		}
		return;
	}
}

template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
AlloyReactionGenerator::addSinks(IndexType i, TTag tag) const
{
	using Species = typename Network::Species;
	using Composition = typename Network::Composition;

	const auto& clReg = this->getCluster(i).getRegion();
	Composition lo = clReg.getOrigin();

	// I
	if (clReg.isSimplex() && lo.isOnAxis(Species::I)) {
		this->addSinkReaction(tag, {i, Network::invalidIndex()});
	}

	// V
	if (clReg.isSimplex() && lo.isOnAxis(Species::V)) {
		this->addSinkReaction(tag, {i, Network::invalidIndex()});
	}
}

template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
AlloyReactionGenerator::addTransforms(IndexType i, IndexType j, TTag tag) const
{
	// Faulted becomes perfect above a size threshold
	auto vSize = this->_clusterData.transitionSizeV();
	auto iSize = this->_clusterData.transitionSizeI();

	using Species = typename Network::Species;
	using Composition = typename Network::Composition;

	const auto& cl1Reg = this->getCluster(i).getRegion();
	Composition lo1 = cl1Reg.getOrigin();

	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo2 = cl2Reg.getOrigin();

	// Interstitial
	if ((lo1[Species::PerfectI] > 0 and lo2[Species::FaultedI] > 0) or
		(lo2[Species::PerfectI] > 0 and lo1[Species::FaultedI] > 0)) {
		// Which is faulted ?
		auto faulted = lo1[Species::FaultedI] > 0 ? i : j;
		auto fReg = lo1[Species::FaultedI] > 0 ? lo1 : lo2;
		auto perfect = lo1[Species::FaultedI] > 0 ? j : i;
		auto pReg = lo1[Species::FaultedI] > 0 ? lo2 : lo1;

		// Need to be the same size
		if (fReg[Species::FaultedI] != pReg[Species::PerfectI])
			return;

		// Large enough
		if (fReg[Species::FaultedI] < iSize)
			return;

		// Add it
		this->addTransformReaction(tag, {faulted, perfect});
	}

	// Vacancy
	if ((lo1[Species::PerfectV] > 0 and lo2[Species::FaultedV] > 0) or
		(lo2[Species::PerfectV] > 0 and lo1[Species::FaultedV] > 0)) {
		// Which is faulted ?
		auto faulted = lo1[Species::FaultedV] > 0 ? i : j;
		auto fReg = lo1[Species::FaultedV] > 0 ? lo1 : lo2;
		auto perfect = lo1[Species::FaultedV] > 0 ? j : i;
		auto pReg = lo1[Species::FaultedV] > 0 ? lo2 : lo1;

		// Need to be the same size
		if (fReg[Species::FaultedV] != pReg[Species::PerfectV])
			return;

		// Large enough
		if (fReg[Species::FaultedV] < vSize)
			return;

		// Add it
		this->addTransformReaction(tag, {faulted, perfect});
	}
}

inline ReactionCollection<AlloyReactionGenerator::Network>
AlloyReactionGenerator::getReactionCollection() const
{
	ReactionCollection<Network> ret(this->_clusterData.gridSize,
		this->getProductionReactions(), this->getDissociationReactions(),
		this->getSinkReactions(), this->getTransformReactions());
	return ret;
}
} // namespace detail

inline detail::AlloyReactionGenerator
AlloyReactionNetwork::getReactionGenerator() const noexcept
{
	return detail::AlloyReactionGenerator{*this};
}
} // namespace network
} // namespace core
} // namespace xolotl
