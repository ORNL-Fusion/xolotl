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
namespace detail
{
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
FeCrReactionGenerator::operator()(IndexType i, IndexType j, TTag tag) const
{
	using Species = typename Network::Species;
	using Composition = typename Network::Composition;
	using AmountType = typename Network::AmountType;

	if (i == j) {
		addSinks(i, tag);
	}

	addTransforms(i, j, tag);

	auto& subpaving = this->getSubpaving();
	auto previousIndex = subpaving.invalidIndex();

	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	// void + void → void
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
		return;
	}

	// void + int → void | int | recombine
	if ((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::I))) {
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
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		auto fSize =
			lo1.isOnAxis(Species::V) ? lo2[Species::Free] : lo1[Species::Free];
		// Compute the product size
		int prodSize = fSize - vSize;
		if (prodSize > 0) {
			//  Find the corresponding cluster
			Composition comp = Composition::zero();
			comp[Species::Free] = prodSize;
			auto fProdId = subpaving.findTileId(comp);
			if (fProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, fProdId});
				// No dissociation
			}
			comp[Species::Free] = 0;
			comp[Species::I] = prodSize;
			auto iProdId = subpaving.findTileId(comp);
			if (iProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, iProdId});
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
			// No product
			this->addProductionReaction(tag, {i, j});
		}
		return;
	}

	// void + loop → loop | complex | trap
	if ((lo1[Species::Loop] > 0 && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2[Species::Loop] > 0)) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		auto lSize =
			lo1.isOnAxis(Species::V) ? lo2[Species::Loop] : lo1[Species::Loop];
		// Compute the product size
		int prodSize = lSize - vSize;
		if (prodSize > 0) {
			//  Find the corresponding cluster
			Composition comp = Composition::zero();
			comp[Species::Loop] = prodSize;
			comp[Species::Trap] = 1;
			auto lProdId = subpaving.findTileId(comp);
			if (lProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, lProdId});
				// No dissociation
			}
			comp[Species::Loop] = 0;
			comp[Species::Complex] = prodSize;
			auto cProdId = subpaving.findTileId(comp);
			if (cProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, cProdId});
				// No dissociation
			}
		}
		else if (prodSize == 0) {
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

	// void + junction → junction | trapped | complex | trap
	if ((lo1[Species::Junction] > 0 && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2[Species::Junction] > 0)) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		auto jSize = lo1.isOnAxis(Species::V) ? lo2[Species::Junction] :
												lo1[Species::Junction];
		// Compute the product size
		int prodSize = jSize - vSize;
		if (prodSize > 0) {
			//  Find the corresponding cluster
			Composition comp = Composition::zero();
			comp[Species::Junction] = prodSize;
			comp[Species::Trap] = 1;
			auto jProdId = subpaving.findTileId(comp);
			if (jProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, jProdId});
				// No dissociation
				return;
			}
			comp[Species::Junction] = 0;
			comp[Species::Trapped] = prodSize;
			auto tProdId = subpaving.findTileId(comp);
			if (tProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, tProdId});
				// No dissociation
				return;
			}
			comp[Species::Trapped] = 0;
			comp[Species::Complex] = prodSize;
			auto cProdId = subpaving.findTileId(comp);
			if (cProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, cProdId});
				// No dissociation
			}
		}
		else if (prodSize == 0) {
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
		// Compute the composition of the new cluster
		auto size = lo1[Species::I] + lo2[Species::I] + lo1[Species::Trapped] +
			lo2[Species::Trapped]; // The other axis should be 0 so it should
								   // work to add up everything
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::Trapped] = size;
		comp[Species::Trap] = 1;
		auto iProdId = subpaving.findTileId(comp);
		if (iProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, iProdId});
			if (lo1[Species::I] == 1 || lo2[Species::I] == 1) {
				this->addDissociationReaction(tag, {iProdId, i, j});
			}
		}
		return;
	}

	// int + free → free
	if ((lo1.isOnAxis(Species::Free) && lo2.isOnAxis(Species::I)) ||
		(lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::Free))) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::I] + lo2[Species::I] + lo1[Species::Free] +
			lo2[Species::Free]; // The other axis should be 0 so it should
								// work to add up everything
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::Free] = size;
		auto pProdId = subpaving.findTileId(comp);
		if (pProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, pProdId});
			if (lo1[Species::I] == 1 || lo2[Species::I] == 1) {
				this->addDissociationReaction(tag, {pProdId, i, j});
			}
		}
		return;
	}

	// int + loop → loop
	if ((lo1[Species::Loop] > 0 && lo2.isOnAxis(Species::I)) ||
		(lo1.isOnAxis(Species::I) && lo2[Species::Loop] > 0)) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::I] + lo2[Species::I] + lo1[Species::Loop] +
			lo2[Species::Loop]; // The other axis should be 0 so it should
								// work to add up everything
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::Loop] = size;
		comp[Species::Trap] = 1;
		auto iProdId = subpaving.findTileId(comp);
		if (iProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, iProdId});
			// No dissociation
		}
		return;
	}

	// int + junction → junction
	if ((lo1[Species::Junction] > 0 && lo2.isOnAxis(Species::I)) ||
		(lo1.isOnAxis(Species::I) && lo2[Species::Junction] > 0)) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::I] + lo2[Species::I] + lo1[Species::Junction] +
			lo2[Species::Junction]; // The other axis should be 0 so it should
									// work to add up everything
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::Junction] = size;
		comp[Species::Trap] = 1;
		auto iProdId = subpaving.findTileId(comp);
		if (iProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, iProdId});
			// No dissociation
		}
		return;
	}

	// free + trapped → junction | trapped
	if ((lo1[Species::Trapped] > 0 && lo2.isOnAxis(Species::Free)) ||
		(lo1.isOnAxis(Species::Free) && lo2[Species::Trapped] > 0)) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::Free] + lo2[Species::Free] +
			lo1[Species::Trapped] +
			lo2[Species::Trapped]; // The other axis should be 0 so it should
								   // work to add up everything

		// Band condition
		double nF = lo1[Species::Free] + lo2[Species::Free];
		double nT = lo1[Species::Trapped] + lo2[Species::Trapped];
		double ratio = std::fabs(nF - nT) / (nF + nT);
		if (ratio < 0.5) {
			// Product is junction
			Composition comp = Composition::zero();
			comp[Species::Junction] = size;
			comp[Species::Trap] = 1;
			auto iProdId = subpaving.findTileId(comp);
			if (iProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, iProdId});
				// No dissociation
			}
		}
		else {
			// Product is trapped
			Composition comp = Composition::zero();
			comp[Species::Trapped] = size;
			comp[Species::Trap] = 1;
			auto iProdId = subpaving.findTileId(comp);
			if (iProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, iProdId});
				// No dissociation
			}
		}
		return;
	}

	// free + trap → trapped
	if ((lo1.isOnAxis(Species::Free) && lo2.isOnAxis(Species::Trap)) ||
		(lo1.isOnAxis(Species::Trap) && lo2.isOnAxis(Species::Free))) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::Free] +
			lo2[Species::Free]; // The other axis should be 0 so it should
								// work to add up everything
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::Trapped] = size;
		comp[Species::Trap] = 1;
		auto tProdId = subpaving.findTileId(comp);
		if (tProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, tProdId});
			this->addDissociationReaction(tag, {tProdId, i, j});
		}

		// Special dissociation
		comp[Species::Trapped] = 0;
		comp[Species::Junction] = size;
		auto jProdId = subpaving.findTileId(comp);
		if (jProdId != subpaving.invalidIndex()) {
			this->addDissociationReaction(tag, {jProdId, i, j});
		}
		return;
	}

	// free + loop → junction | trapped | loop
	if ((lo1[Species::Loop] > 0 && lo2.isOnAxis(Species::Free)) ||
		(lo1.isOnAxis(Species::Free) && lo2[Species::Loop] > 0)) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::Free] + lo2[Species::Free] +
			lo1[Species::Loop] +
			lo2[Species::Loop]; // The other axis should be 0 so it should
								// work to add up everything

		// Band condition
		double nF = lo1[Species::Free] + lo2[Species::Free];
		double nL = lo1[Species::Loop] + lo2[Species::Loop];
		double ratio = std::fabs(nF - nL) / (nF + nL);
		if (ratio < 0.5) {
			// Product is junction
			Composition comp = Composition::zero();
			comp[Species::Junction] = size;
			comp[Species::Trap] = 1;
			auto iProdId = subpaving.findTileId(comp);
			if (iProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, iProdId});
				// No dissociation
			}
		}
		else {
			// Restrictions
			if (nF < nL) {
				// Product is loop
				Composition comp = Composition::zero();
				comp[Species::Loop] = size;
				comp[Species::Trap] = 1;
				auto iProdId = subpaving.findTileId(comp);
				if (iProdId != subpaving.invalidIndex()) {
					this->addProductionReaction(tag, {i, j, iProdId});
					// No dissociation
				}
			}
			else {
				// Product is trapped
				Composition comp = Composition::zero();
				comp[Species::Trapped] = size;
				comp[Species::Trap] = 1;
				auto iProdId = subpaving.findTileId(comp);
				if (iProdId != subpaving.invalidIndex()) {
					this->addProductionReaction(tag, {i, j, iProdId});
					// No dissociation
				}
			}
		}
		return;
	}

	// free + junction → junction
	if ((lo1[Species::Junction] > 0 && lo2.isOnAxis(Species::Free)) ||
		(lo1.isOnAxis(Species::Free) && lo2[Species::Junction] > 0)) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::Free] + lo2[Species::Free] +
			lo1[Species::Junction] +
			lo2[Species::Junction]; // The other axis should be 0 so it should
									// work to add up everything

		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::Junction] = size;
		comp[Species::Trap] = 1;
		auto jProdId = subpaving.findTileId(comp);
		if (jProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, jProdId});
			// No dissociation
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
	if (clReg.isSimplex() && lo.isOnAxis(Species::Free)) {
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

	const auto& clReg = this->getCluster(i).getRegion();
	Composition lo = clReg.getOrigin();

	const auto& prReg = this->getCluster(j).getRegion();
	Composition loPr = prReg.getOrigin();

	// The first reactant need to be Junction
	if (clReg.isSimplex() && lo[Species::Junction] == 0)
		return;

	// The second on should be Trapped or Loop
	if (prReg.isSimplex() && loPr[Species::Trapped] == 0 &&
		loPr[Species::Loop] == 0)
		return;

	// They need to be the same size
	if (loPr[Species::Trapped] == lo[Species::Junction] ||
		loPr[Species::Loop] == lo[Species::Junction]) {
		this->addTransformReaction(tag, {i, j});

		//		std::cout << lo[static_cast<int>(Species::V)] << "-" <<
		//lo[static_cast<int>(Species::I)] << "-" <<
		//lo[static_cast<int>(Species::Free)] << "-"
		//				 << lo[static_cast<int>(Species::Trapped)] << "-" <<
		//lo[static_cast<int>(Species::Junction)] << "-" <<
		//lo[static_cast<int>(Species::Complex)] << "-"
		//				 << lo[static_cast<int>(Species::Loop)] << " + "  <<
		//loPr[static_cast<int>(Species::V)] << "-" <<
		//loPr[static_cast<int>(Species::I)] << "-" <<
		//loPr[static_cast<int>(Species::Free)] << "-"
		//				 << loPr[static_cast<int>(Species::Trapped)] << "-" <<
		//loPr[static_cast<int>(Species::Junction)] << "-" <<
		//loPr[static_cast<int>(Species::Complex)] << "-"
		//				 << loPr[static_cast<int>(Species::Loop)] << std::endl;
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
