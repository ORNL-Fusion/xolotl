#pragma once

#include <xolotl/core/network/detail/impl/ConstantReactionGenerator.tpp>
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
	if (not map["SSBM"])
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
	if (this->_enableSSBM) {
		auto clusterDataMirror = this->getClusterDataMirror();
		
		// Get the concentrations on the host
		auto dConcs = Kokkos::subview(concentrations,
			std::make_pair(
				clusterDataMirror.bubbleId(), clusterDataMirror.bubbleId() + 2));
		auto heConcs = create_mirror_view(dConcs);
		
		deep_copy(heConcs, dConcs);
	
		// Compute the average composition of each defect
			auto conc = heConcs(0);
			auto avComp = heConcs(1) / conc;
			if (conc == 0.0 )
				avComp = 0.0;
		// Void
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
	if (this->_enableSSBM) {
		auto clusterDataMirror = this->getClusterDataMirror();

		// Get the concentrations on the host SSBM
		auto dConcs = Kokkos::subview(concentrations,
			std::make_pair(
				clusterDataMirror.bubbleId(),clusterDataMirror.bubbleId() + 2));
		auto heConcs = create_mirror_view(dConcs);
		deep_copy(heConcs, dConcs);

		// Compute the average composition for V
			auto conc = heConcs(0);
			auto avComp = 1; //heConcs(1) / conc;
		//	if (conc <= 0.0 || heConcs(1) <= 0)
		//		avComp = 0.0;
		// Void
			this->_clusterData.h_view().setBubbleAvRad(util::max(0.0,
				computeBubbleRadius(
					avComp, clusterDataMirror.latticeParameter())));
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
	
	constexpr auto speciesRange = NetworkType::getSpeciesRange();
	constexpr auto speciesNoI = NetworkType::getSpeciesRangeNoI();

	auto& subpaving = this->getSubpaving();
	auto previousIndex = subpaving.invalidIndex();

	// Sinks
	if (i == j) {
		addSinks(i, tag);
	}

	// Add the traps between He and Traps
	// addTraps(i, j, tag);
	
	auto numClusters = this->getNumberOfClusters();

	// SSBM

	// Large bubble reactions
	if (this->_clusterData.enableSSBM())
		addSingleSizeReactions(i, j, tag);
	
	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();


	// V + V = V
	if (lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::V) &&
		diffusionFactor(i) > 0.0 && diffusionFactor(j) > 0.0) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::V] + lo2[Species::V];
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::V] = size;
		auto vProdId = subpaving.findTileId(comp);
		if (vProdId != subpaving.invalidIndex() && vProdId != previousIndex) {
			this->addProductionReaction(tag, {i, j, vProdId});
			if (lo1[Species::V] == 1 || lo2[Species::V] == 1) {
				this->addDissociationReaction(tag, {vProdId, i, j});
			}
			previousIndex = vProdId;
		}
		return;
	}

	// vac + int = vac | int | recombine
	if (((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
			(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::I)))) {
		// Find out which one is which
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];
		// Vac can be grouped
		auto vReg = lo1.isOnAxis(Species::V) ? cl1Reg : cl2Reg;
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
				if (iProdId != subpaving.invalidIndex() &&
					iProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, iProdId});
					previousIndex = iProdId;
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

	// He + V = HeV
	if ((lo1.isOnAxis(Species::He) and lo2.isOnAxis(Species::V)) or
		(lo1.isOnAxis(Species::V) and lo2.isOnAxis(Species::He))) {
		// Void can be grouped
		auto minSize = lo1[Species::V] + lo2[Species::V];
		auto maxSize = hi1[Species::V] + hi2[Species::V] - 2;
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::He] = lo1[Species::He] + lo2[Species::He];
			comp[Species::V] = k;
			auto vProdId = subpaving.findTileId(comp);
			if (vProdId != subpaving.invalidIndex() &&
				vProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, vProdId});

				// No dissociation
				previousIndex = vProdId;
			}
		}
		return;
	}

	// He + HeV = HeV
	if (lo1[Species::He] > 0 and lo2[Species::He] > 0) {
		// HeV can be grouped
		auto minVSize = lo1[Species::V] + lo2[Species::V];
		auto maxVSize = hi1[Species::V] + hi2[Species::V] - 2;
		auto minHeSize = lo1[Species::He] + lo2[Species::He];
		auto maxHeSize = hi1[Species::He] + hi2[Species::He] - 2;
		// Find the corresponding clusters
		for (auto k = minHeSize; k <= maxHeSize; k++) {
			for (auto l = minVSize; l <= maxVSize; l++) {
				Composition comp = Composition::zero();
				comp[Species::He] = k;
				comp[Species::V] = l;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() &&
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					// No dissociation
					previousIndex = vProdId;
				}
				if (vProdId == subpaving.invalidIndex()) {
					// Trap mutation
					comp[Species::V] = l + 1;
					vProdId = subpaving.findTileId(comp);
					if (vProdId != subpaving.invalidIndex() &&
						vProdId != previousIndex) {
						// Get the I product
						Composition iComp = Composition::zero();
						iComp[Species::I] = 1;
						auto iProdId = subpaving.findTileId(iComp);
						if (iProdId != subpaving.invalidIndex()) {
							this->addProductionReaction(
								tag, {i, j, vProdId, iProdId});
						}
					}
				}
			}
		}
		return;
	}

	// V + HeV = HeV
	if (lo1.isOnAxis(Species::V) or lo2.isOnAxis(Species::V)) {
		// HeV can be grouped
		auto minVSize = lo1[Species::V] + lo2[Species::V];
		auto maxVSize = hi1[Species::V] + hi2[Species::V] - 2;
		auto minHeSize = lo1[Species::He] + lo2[Species::He];
		auto maxHeSize = hi1[Species::He] + hi2[Species::He] - 2;
		// Find the corresponding clusters
		for (auto l = minVSize; l <= maxVSize; l++) {
			for (auto k = minHeSize; k <= maxHeSize; k++) {
				Composition comp = Composition::zero();
				comp[Species::He] = k;
				comp[Species::V] = l;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() &&
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					// No dissociation
					previousIndex = vProdId;
				}
			}
		}
		return;
	}


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
	const auto& clReg = this->getCluster(i).getRegion();
        Composition lo = clReg.getOrigin();	

	// Check reaction with largest bubble
	if (clReg.isSimplex()){
	if (i == j) {

		// V case
		if (lo.isOnAxis(Species::V)) {
			// V_k + B -> B
			this->addProductionReaction(tag, {i, bubbleId, bubbleId});
		}
		// I case
		else if (lo.isOnAxis(Species::I)) {
			// I_k + B -> B
			this->addProductionReaction(tag, {i, bubbleId, bubbleId});
		}

		// He case
		else if (lo.isOnAxis(Species::He)) {
			// He_k + B -> B
			this->addProductionReaction(tag, {i, bubbleId, bubbleId});
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
	auto largestHeSize = hiLargest[Species::He] - 1;

	// He_a + He_bV -> B
	if (hi1[Species::He] + hi2[Species::He] - 2 > largestHeSize) {
		this->addProductionReaction(tag, {i, j, bubbleId});
	}

	// V_a + HeV_b -> B
	if (hi1[Species::V] + hi2[Species::V] - 2 > largestVSize) {
		this->addProductionReaction(tag, {i, j, bubbleId});

	}


	// I_a + B -> HeV_b
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
