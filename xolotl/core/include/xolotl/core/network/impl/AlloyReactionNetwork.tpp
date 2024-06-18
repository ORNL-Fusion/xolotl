#pragma once

#include <xolotl/core/network/detail/impl/ConstantReactionGenerator.tpp>
#include <xolotl/core/network/detail/impl/SinkReactionGenerator.tpp>
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
	// Get the diffusion factors
	auto diffusionFactor = this->_clusterData.diffusionFactor;

	using Species = typename Network::Species;
	using Composition = typename Network::Composition;
	using AmountType = typename Network::AmountType;

	if (i == j) {
		if (diffusionFactor(i) != 0.0)
			addSinks(i, tag);

		if (this->_constantConnsRows.extent(0) > 0) {
			// Look for the entry
			for (auto k = this->_constantConnsRows(i);
				 k < this->_constantConnsRows(i + 1); k++) {
				if (this->_constantConnsEntries(k) == this->_numDOFs) {
					this->addConstantReaction(
						tag, {i, Network::invalidIndex()});
					break;
				}
			}
		}
	}

	// Add every possibility
	if (this->_constantConnsRows.extent(0) > 0) {
		// Look for the entry
		for (auto k = this->_constantConnsRows(i);
			 k < this->_constantConnsRows(i + 1); k++) {
			if (this->_constantConnsEntries(k) == j) {
				this->addConstantReaction(tag, {i, j});
				break;
			}
		}
	}
	if (j != i) {
		if (this->_constantConnsRows.extent(0) > 0) {
			// Look for the entry
			for (auto k = this->_constantConnsRows(j);
				 k < this->_constantConnsRows(j + 1); k++) {
				if (this->_constantConnsEntries(k) == i) {
					this->addConstantReaction(tag, {j, i});
					break;
				}
			}
		}
	}

	auto& subpaving = this->getSubpaving();
	auto previousIndex = subpaving.invalidIndex();

	// Check the diffusion factor
	if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
		return;
	}

	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	// vac + vac = vac | void and faulted
	if (lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::V)) {
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
		comp[Species::V] = 0;
		comp[Species::Void] = size;
		auto voidProdId = subpaving.findTileId(comp);
		if (voidProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, voidProdId});
			if (lo1[Species::V] == 1 || lo2[Species::V] == 1) {
				this->addDissociationReaction(tag, {voidProdId, i, j});
			}
		}
		comp[Species::Void] = 0;
		comp[Species::Faulted] = size;
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
	if ((lo1.isOnAxis(Species::Void) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::Void))) {
		// Void can be grouped
		auto minSize = lo1[Species::V] + lo2[Species::V] + lo1[Species::Void] +
			lo2[Species::Void];
		auto maxSize = lo1[Species::V] + lo2[Species::V] + hi1[Species::Void] +
			hi2[Species::Void] - 2; // The other axis should be 0 so it should
									// work to add up everything
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::Void] = k;
			auto voidProdId = subpaving.findTileId(comp);
			if (voidProdId != subpaving.invalidIndex() &&
				voidProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, voidProdId});
				if (lo1[Species::V] == 1 || lo2[Species::V] == 1) {
					this->addDissociationReaction(tag, {voidProdId, i, j});
				}
				previousIndex = voidProdId;
			}
		}
		return;
	}

	// vac + faulted = faulted
	if ((lo1.isOnAxis(Species::Faulted) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::Faulted))) {
		// Faulted can be grouped
		auto minSize = lo1[Species::V] + lo2[Species::V] +
			lo1[Species::Faulted] + lo2[Species::Faulted];
		auto maxSize = lo1[Species::V] + lo2[Species::V] +
			hi1[Species::Faulted] + hi2[Species::Faulted] -
			2; // The other axis should be 0 so it should work
			   // to add up everything
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::Faulted] = k;
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

	// vac + perfect = perfect | int | recombine | vac
	if (((lo1.isOnAxis(Species::Perfect) && lo2.isOnAxis(Species::V)) ||
			(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::Perfect)))) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		// Perfect can be grouped
		auto pReg = lo1.isOnAxis(Species::Perfect) ? cl1Reg : cl2Reg;
		for (auto k : makeIntervalRange(pReg[Species::Perfect])) {
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
				comp[Species::Perfect] = prodSize;
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

	// vac + frank = frank | int | recombine | vac
	if (((lo1.isOnAxis(Species::Frank) && lo2.isOnAxis(Species::V)) ||
			(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::Frank)))) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		// Frank can be grouped
		auto fReg = lo1.isOnAxis(Species::Frank) ? cl1Reg : cl2Reg;
		for (auto k : makeIntervalRange(fReg[Species::Frank])) {
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
				comp[Species::Frank] = prodSize;
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

	// int + int = int | frank and perfect
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
		comp[Species::Frank] = size;
		auto fProdId = subpaving.findTileId(comp);
		if (fProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, fProdId});
			// No dissociation
		}
		comp[Species::Frank] = 0;
		comp[Species::Perfect] = size;
		auto pProdId = subpaving.findTileId(comp);
		if (pProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, pProdId});
			// No dissociation
		}
		return;
	}

	// int + frank = frank
	if ((lo1.isOnAxis(Species::Frank) && lo2.isOnAxis(Species::I)) ||
		(lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::Frank))) {
		// Frank can be grouped
		auto minSize = lo1[Species::I] + lo2[Species::I] + lo1[Species::Frank] +
			lo2[Species::Frank];
		auto maxSize = lo1[Species::I] + lo2[Species::I] + hi1[Species::Frank] +
			hi2[Species::Frank] - 2; // The other axis should be 0 so it should
									 // work to add up everything
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::Frank] = k;
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

	// int + perfect = perfect
	if ((lo1.isOnAxis(Species::Perfect) && lo2.isOnAxis(Species::I)) ||
		(lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::Perfect))) {
		// Perfect can be grouped
		auto minSize = lo1[Species::I] + lo2[Species::I] +
			lo1[Species::Perfect] + lo2[Species::Perfect];
		auto maxSize = lo1[Species::I] + lo2[Species::I] +
			hi1[Species::Perfect] + hi2[Species::Perfect] -
			2; // The other axis should be 0 so it should
			   // work to add up everything
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::Perfect] = k;
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

	// int + faulted = faulted | vac
	if (((lo1.isOnAxis(Species::Faulted) && lo2.isOnAxis(Species::I)) ||
			(lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::Faulted)))) {
		// Find out which one is which
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];
		// Faulted can be grouped
		auto fReg = lo1.isOnAxis(Species::Faulted) ? cl1Reg : cl2Reg;
		for (auto k : makeIntervalRange(fReg[Species::Faulted])) {
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
				comp[Species::Faulted] = prodSize;
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

	// int + void = void | vac
	if (((lo1.isOnAxis(Species::Void) && lo2.isOnAxis(Species::I)) ||
			(lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::Void)))) {
		// Find out which one is which
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];
		// Void can be grouped
		auto vReg = lo1.isOnAxis(Species::Void) ? cl1Reg : cl2Reg;
		for (auto k : makeIntervalRange(vReg[Species::Void])) {
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
				comp[Species::Void] = prodSize;
				auto fProdId = subpaving.findTileId(comp);
				if (fProdId != subpaving.invalidIndex() &&
					fProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, fProdId});
					//					if (lo1[Species::I] == 1 ||
					// lo2[Species::I]
					//== 1) {
					// this->addDissociationReaction(tag, {fProdId, i, j});
					//					}
					previousIndex = fProdId;
				}
			}
		}
		return;
	}

	// perfect + faulted = faulted | vac | int | perfect | recombine
	if (((lo1.isOnAxis(Species::Perfect) && lo2.isOnAxis(Species::Faulted)) ||
			(lo1.isOnAxis(Species::Faulted) &&
				lo2.isOnAxis(Species::Perfect)))) {
		// Perfect and Faulted can be grouped
		auto minSize = lo1[Species::Faulted] + lo2[Species::Faulted] -
			hi1[Species::Perfect] - hi2[Species::Perfect] + 2;
		auto maxSize = hi1[Species::Faulted] + hi2[Species::Faulted] -
			lo1[Species::Perfect] - lo2[Species::Perfect] -
			2; // The other axis should be 0 so it should
			   // work to add up everything

		for (auto k = minSize; k < maxSize; k++) {
			// 3 cases
			if (k > 0) {
				// Looking for V cluster
				Composition comp = Composition::zero();
				comp[Species::V] = k;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() &&
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					previousIndex = vProdId;
					// No dissociation
				}
				comp[Species::V] = 0;
				comp[Species::Faulted] = k;
				auto fProdId = subpaving.findTileId(comp);
				if (fProdId != subpaving.invalidIndex() &&
					fProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, fProdId});
					previousIndex = fProdId;
					// No dissociation
				}
			}
			else if (k < 0) {
				// Looking for I cluster
				Composition comp = Composition::zero();
				comp[Species::I] = -k;
				auto iProdId = subpaving.findTileId(comp);
				if (iProdId != subpaving.invalidIndex() &&
					iProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, iProdId});
					previousIndex = iProdId;
					// No dissociation
				}
				comp[Species::I] = 0;
				comp[Species::Perfect] = -k;
				auto pProdId = subpaving.findTileId(comp);
				if (pProdId != subpaving.invalidIndex() &&
					pProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, pProdId});
					previousIndex = pProdId;
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

	// perfect + frank = frank
	if ((lo1.isOnAxis(Species::Frank) && lo2.isOnAxis(Species::Perfect)) ||
		(lo1.isOnAxis(Species::Perfect) && lo2.isOnAxis(Species::Frank))) {
		// Perfect and Frank can be grouped
		auto minSize = lo1[Species::Perfect] + lo2[Species::Perfect] +
			lo1[Species::Frank] + lo2[Species::Frank];
		auto maxSize = hi1[Species::Perfect] + hi2[Species::Perfect] +
			hi1[Species::Frank] + hi2[Species::Frank] -
			4; // The other axis should be 0 so it should work
			   // to add up everything
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::Frank] = k;
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

	// perfect + void = void | vac | int | perfect | recombine
	if (((lo1.isOnAxis(Species::Void) && lo2.isOnAxis(Species::Perfect)) ||
			(lo1.isOnAxis(Species::Perfect) && lo2.isOnAxis(Species::Void)))) {
		// Perfect and Faulted can be grouped
		auto minSize = lo1[Species::Void] + lo2[Species::Void] -
			hi1[Species::Perfect] - hi2[Species::Perfect] + 2;
		auto maxSize = hi1[Species::Void] + hi2[Species::Void] -
			lo1[Species::Perfect] - lo2[Species::Perfect] -
			2; // The other axis should be 0 so it should
			   // work to add up everything

		for (auto k = minSize; k < maxSize; k++) {
			// 3 cases
			if (k > 0) {
				// Looking for V cluster
				Composition comp = Composition::zero();
				comp[Species::V] = k;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() &&
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					previousIndex = vProdId;
					// No dissociation
				}
				comp[Species::V] = 0;
				comp[Species::Void] = k;
				vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() &&
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					previousIndex = vProdId;
					// No dissociation
				}
			}
			else if (k < 0) {
				// Looking for I cluster
				Composition comp = Composition::zero();
				comp[Species::I] = -k;
				auto iProdId = subpaving.findTileId(comp);
				if (iProdId != subpaving.invalidIndex() &&
					iProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, iProdId});
					previousIndex = iProdId;
					// No dissociation
				}
				comp[Species::I] = 0;
				comp[Species::Perfect] = -k;
				auto pProdId = subpaving.findTileId(comp);
				if (pProdId != subpaving.invalidIndex() &&
					pProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, pProdId});
					previousIndex = pProdId;
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

	// perfect + perfect = perfect
	if (lo1.isOnAxis(Species::Perfect) && lo2.isOnAxis(Species::Perfect)) {
		// Perfect can be grouped
		auto minSize = lo1[Species::Perfect] + lo2[Species::Perfect];
		auto maxSize = hi1[Species::Perfect] + hi2[Species::Perfect] -
			2; // The other axis should be 0 so it should work
			   // to add up everything
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::Perfect] = k;
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

	// Perfect
	if (clReg.isSimplex() && lo.isOnAxis(Species::Perfect)) {
		this->addSinkReaction(tag, {i, Network::invalidIndex()});
	}
}

inline ReactionCollection<AlloyReactionGenerator::Network>
AlloyReactionGenerator::getReactionCollection() const
{
	ReactionCollection<Network> ret(this->_clusterData.gridSize,
		this->getProductionReactions(), this->getDissociationReactions(),
		this->getSinkReactions(), this->getConstantReactions());
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
