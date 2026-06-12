#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/detail/impl/SinkReactionGenerator.tpp>
#include <xolotl/core/network/detail/impl/TrapReactionGenerator.tpp>
#include <xolotl/core/network/impl/FeClusterGenerator.tpp>
#include <xolotl/core/network/impl/FeReaction.tpp>
#include <xolotl/core/network/impl/ReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{

void
FeReactionNetwork::initializeExtraDOFs(const options::IOptions& options)
{
	auto map = options.getProcesses();
	if (not map["largeBubble"])
		return;

	largestClusterId = checkLargestClusterId();

	this->_clusterData.h_view().setBubbleId(this->_numDOFs);
	this->_clusterData.h_view().setVoidAvId(this->_numDOFs + 1);
	this->_clusterData.h_view().setHeAvId(this->_numDOFs + 2);
	this->_numDOFs += 3;
}

void
FeReactionNetwork::computeFluxesPreProcess(ConcentrationsView concentrations,
	FluxesView fluxes, IndexType gridIndex, double surfaceDepth, double spacing)
{
	if (this->_enableLargeBubble) {
		auto clusterDataMirror = this->getClusterDataMirror();

		// Get the concentrations on the host
		auto dConcs = Kokkos::subview(concentrations,
			std::make_pair(
				clusterDataMirror.bubbleId(), clusterDataMirror.bubbleId() + 2));
		auto heConcs = create_mirror_view(dConcs);
		deep_copy(heConcs, dConcs);

		// Compute the average composition for V

			auto conc = heConcs(0);
			auto avComp = heConcs(1) / conc;
			if (conc == 0.0)
				avComp = 0.0;
			// Compute and save the radius from that
				this->_clusterData.h_view().setBubbleAvRad(util::max(0.0,
					computeBubbleRadius(
						avComp, clusterDataMirror.latticeParameter())));
	}
}

void
FeReactionNetwork::computePartialsPreProcess(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex, double surfaceDepth, double spacing)
{
	if (this->_enableLargeBubble) {
		auto clusterDataMirror = this->getClusterDataMirror();

		// Get the concentrations on the host SSBM
		auto dConcs = Kokkos::subview(concentrations,
			std::make_pair(
				clusterDataMirror.bubbleId(),clusterDataMirror.bubbleId() + 2));
		auto heConcs = create_mirror_view(dConcs);
		deep_copy(heConcs, dConcs);

		// Compute the average composition for V
			auto conc = heConcs(0);
			auto avComp = heConcs(1) / conc;
			if (conc == 0.0)
				avComp = 0.0;
			// Void
				this->_clusterData.h_view().setBubbleAvRad(util::max(0.0,
					computeBubbleRadius(avComp, clusterDataMirror.latticeParameter())));
		}
}


namespace detail
{
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
FeReactionGenerator::operator()(IndexType i, IndexType j, TTag tag) const
{
		// Check the diffusion factors
	auto diffusionFactor = this->_clusterData.diffusionFactor;
	if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
		return;
	}

	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;

	constexpr auto species = NetworkType::getSpeciesRange();
	constexpr auto speciesNoI = NetworkType::getSpeciesRangeNoI();

	// Add the sinks
	if (i == j) {
		addSinks(i, tag);
	}

	// SSBM
	if (this->_clusterData.enableSSBM())
		addSingleSizeReactions(i, j, tag);
	
	// Add the traps between He and Traps


	auto numClusters = this->getNumberOfClusters();

	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	auto& subpaving = this->getSubpaving();
	auto previousIndex = subpaving.invalidIndex();

	// Special case for I + I
	if (cl1Reg.isSimplex() && cl2Reg.isSimplex() && lo1.isOnAxis(Species::I) &&
		lo2.isOnAxis(Species::I)) {
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
		return;
	}

	// Special case for I + V
	if ((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::I))) {
		// Compute the largest possible product and the smallest one
		int largestProd = (int)hi1[Species::V] + (int)hi2[Species::V] - 2 -
			(int)lo1[Species::I] - (int)lo2[Species::I];
		int smallestProd = (int)lo1[Species::V] + (int)lo2[Species::V] -
			(int)hi1[Species::I] - (int)hi2[Species::I] + 2;
		// Loop on the products
		for (int prodSize = smallestProd; prodSize <= largestProd; prodSize++) {
			// 3 cases
			if (prodSize > 0) {
				// Looking for V cluster
				Composition comp = Composition::zero();
				comp[Species::V] = prodSize;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() &&
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					// No dissociation
					previousIndex = vProdId;
				}
			}
			else if (prodSize < 0) {
				// Looking for I cluster
				Composition comp = Composition::zero();
				comp[Species::I] = -prodSize;
				auto iProdId = subpaving.findTileId(comp);
				if (iProdId != subpaving.invalidIndex() &&
					iProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, iProdId});
					// No dissociation
					previousIndex = iProdId;
				}
			}
			else {
				// No product
				this->addProductionReaction(tag, {i, j});
			}
		}
		return;
	}

	// General case
	constexpr auto numSpeciesNoI = NetworkType::getNumberOfSpeciesNoI();
	using BoundsArray =
		Kokkos::Array<Kokkos::pair<AmountType, AmountType>, numSpeciesNoI>;
	plsm::EnumIndexed<BoundsArray, Species> bounds;
	// Loop on the species
	for (auto l : species) {
		auto low = lo1[l] + lo2[l];
		auto high = hi1[l] + hi2[l] - 2;
		// Special case for I
		if (l == Species::I) {
			bounds[Species::V].first -= high;
			bounds[Species::V].second -= low;
		}
		else {
			bounds[l] = {low, high};
		}
	}

	// Look for potential product
	IndexType nProd = 0;
	for (IndexType k = 0; k < numClusters; ++k) {
		// Get the composition
		const auto& prodReg = this->getCluster(k).getRegion();
		bool isGood = true;
		// Loop on the species
		for (auto l : speciesNoI) {
			if (prodReg[l()].begin() > bounds[l()].second) {
				isGood = false;
				break;
			}
			if (prodReg[l()].end() - 1 < bounds[l()].first) {
				isGood = false;
				break;
			}
		}

		if (isGood) {
			// Increase nProd
			nProd++;
			this->addProductionReaction(tag, {i, j, k});

			// Loop on the species
			bool isOnAxis1 = false, isOnAxis2 = false;
			for (auto l : species) {
				if (lo1.isOnAxis(l()) && lo1[l()] == 1)
					isOnAxis1 = true;
				if (lo2.isOnAxis(l()) && lo2[l()] == 1)
					isOnAxis2 = true;
			}
			if (isOnAxis1 || isOnAxis2) {
				if (lo1.isOnAxis(Species::He) && lo2.isOnAxis(Species::He) &&
					lo1[Species::He] == 1 && lo2[Species::He] == 1) {
					continue;
				}

				this->addDissociationReaction(tag, {k, i, j});
			}
		}
	}
}

template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
FeReactionGenerator::addSinks(IndexType i, TTag tag) const
{
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	const auto& clReg = this->getCluster(i).getRegion();
	Composition lo = clReg.getOrigin();

	// I
	if (clReg.isSimplex() && lo.isOnAxis(Species::I)) {
		this->addSinkReaction(tag, {i, NetworkType::invalidIndex()});
	}

	// V
	if (clReg.isSimplex() && lo.isOnAxis(Species::V) && lo[Species::V] < 5) {
		this->addSinkReaction(tag, {i, NetworkType::invalidIndex()});
	}

	// He
	if (clReg.isSimplex() && lo.isOnAxis(Species::He) && lo[Species::He] < 4) {
		this->addSinkReaction(tag, {i, NetworkType::invalidIndex()});
	}
}

template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
FeReactionGenerator::addSingleSizeReactions(
	IndexType i, IndexType j, TTag tag) const
{
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	IndexType bubbleId = this->_clusterData.bubbleId();

	if (i == j) {
		const auto& clReg = this->getCluster(i).getRegion();
		Composition lo = clReg.getOrigin();

		// Check reaction with largest bubble
		if (not clReg.isSimplex())
			return;

		// V case
		if (lo.isOnAxis(Species::V)) {
			// V_k + B -> B
			this->addProductionReaction(tag, {i, bubbleId, bubbleId});
		}
		// I case
		if (lo.isOnAxis(Species::I)) {
			// I_k + B -> B
			this->addProductionReaction(tag, {i, bubbleId, bubbleId});
		}

		// H case
		
			if (lo.isOnAxis(Species::He)) {
				// H_k + B -> B

				// Only add trap mutation so that at run time it adds the I
				// concentration if needed.
				auto& subpaving = this->getSubpaving();
				Composition comp = Composition::zero();
				comp[Species::I] = 1;
				auto iClusterId = subpaving.findTileId(comp);
				if (iClusterId == NetworkType::invalidIndex()) {
					this->addProductionReaction(tag, {i, bubbleId, bubbleId});
				}
				else {
					this->addProductionReaction(
						tag, {i, bubbleId, bubbleId, iClusterId});
				}
			}
		
	}
	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	// Find the edge of the phase space
	const auto& largestReg = this->getCluster(largestClusterId).getRegion();
	Composition hiLargest = largestReg.getUpperLimitPoint();
	auto largestVSize = hiLargest[Species::V] - 1;
	auto largestImpSize = hiLargest[Species::He] - 1;

	// H_a + H_bV -> B
		if (hi1[Species::He] + hi2[Species::He] - 2 > largestImpSize) {
			this->addProductionReaction(tag, {i, j, bubbleId});
		}

	// V_a + HV_b -> B
	if (hi1[Species::V] + hi2[Species::V] - 2 > largestVSize) {
		this->addProductionReaction(tag, {i, j, bubbleId});
	}

	// I_a + B -> HV_b
	if ((lo1.isOnAxis(Species::I) and lo2[Species::V] > 0) or
		(lo1[Species::V] > 0 and lo2.isOnAxis(Species::I))) {
		// It should be around the largest size value
		if (hi1[Species::V] + hi2[Species::V] + hi1[Species::I] +
				hi2[Species::I] - 4 >
			largestVSize) {
			// Need to know which one is I
			auto iId = lo1[Species::I] > 0 ? i : j;
			auto vId = lo1[Species::I] > 0 ? j : i;
			this->addProductionReaction(tag, {iId, bubbleId, vId});
		}
	}
}


template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
FeReactionGenerator::addTraps(IndexType i, IndexType j, TTag tag) const
{
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Look for He 1 and trap
	const auto& clReg1 = this->getCluster(i).getRegion();
	Composition lo1 = clReg1.getOrigin();
	const auto& clReg2 = this->getCluster(j).getRegion();
	Composition lo2 = clReg2.getOrigin();

	if (clReg1.isSimplex() and clReg2.isSimplex()) {
		if ((lo1.isOnAxis(Species::He) and lo2.isOnAxis(Species::Trap)) or
			(lo2.isOnAxis(Species::He) and lo1.isOnAxis(Species::Trap))) {
			// Which one is which?
			auto trapId = lo1.isOnAxis(Species::Trap) ? i : j;
			auto heId = lo1.isOnAxis(Species::Trap) ? j : i;

			// Only single helium can trap
			Composition loHe = this->getCluster(heId).getRegion().getOrigin();
			if (loHe[Species::He] == 1) {
				this->addTrapReaction(tag, {heId, trapId});
			}
		}
	}
}

inline ReactionCollection<FeReactionGenerator::NetworkType>
FeReactionGenerator::getReactionCollection() const
{
	ReactionCollection<NetworkType> ret(this->_clusterData.gridSize,
		this->_clusterData.numClusters, this->_enableReadRates,
		this->getProductionReactions(), this->getDissociationReactions(),
		this->getSinkReactions(), this->getTrapReactions());
	return ret;
}
} // namespace detail

inline detail::FeReactionGenerator
FeReactionNetwork::getReactionGenerator() const noexcept
{
	return detail::FeReactionGenerator{*this};
}
} // namespace network
} // namespace core
} // namespace xolotl
