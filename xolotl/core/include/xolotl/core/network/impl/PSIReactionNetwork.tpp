#pragma once

#include <set>

#include <xolotl/core/network/detail/PSITrapMutation.h>
#include <xolotl/core/network/detail/impl/BurstingReactionGenerator.tpp>
#include <xolotl/core/network/detail/impl/SinkReactionGenerator.tpp>
#include <xolotl/core/network/detail/impl/TrapMutationClusterData.tpp>
#include <xolotl/core/network/detail/impl/TrapMutationReactionGenerator.tpp>
#include <xolotl/core/network/impl/PSIClusterGenerator.tpp>
#include <xolotl/core/network/impl/PSIReaction.tpp>
#include <xolotl/core/network/impl/ReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TSpeciesEnum>
PSIReactionNetwork<TSpeciesEnum>::PSIReactionNetwork(const Subpaving& subpaving,
	IndexType gridSize, const options::IOptions& options) :
	Superclass(subpaving, gridSize, options),
	_tmHandler(psi::getTrapMutationHandler(
		this->_enableTrapMutation, options.getMaterial()))
{
}

template <typename TSpeciesEnum>
PSIReactionNetwork<TSpeciesEnum>::PSIReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts,
	const std::vector<SubdivisionRatio>& subdivisionRatios, IndexType gridSize,
	const options::IOptions& options) :
	Superclass(maxSpeciesAmounts, subdivisionRatios, gridSize, options),
	_tmHandler(psi::getTrapMutationHandler(
		this->_enableTrapMutation, options.getMaterial()))
{
}

template <typename TSpeciesEnum>
PSIReactionNetwork<TSpeciesEnum>::PSIReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts, IndexType gridSize,
	const options::IOptions& options) :
	Superclass(maxSpeciesAmounts, gridSize, options),
	_tmHandler(psi::getTrapMutationHandler(
		this->_enableTrapMutation, options.getMaterial()))
{
}

template <typename TSpeciesEnum>
void
PSIReactionNetwork<TSpeciesEnum>::initializeExtraClusterData(
	const options::IOptions& options)
{
	if (!this->_enableTrapMutation) {
		return;
	}

	this->_clusterData.h_view().extraData.trapMutationData.initialize();
	this->copyClusterDataView();
	this->invalidateDataMirror();
}

template <typename TSpeciesEnum>
void
PSIReactionNetwork<TSpeciesEnum>::initializeExtraDOFs(
	const options::IOptions& options)
{
	auto map = options.getProcesses();
	if (not map["largeBubble"])
		return;

	largestClusterId = checkLargestClusterId();

	this->_clusterData.h_view().setBubbleId(this->_numDOFs);
	this->_clusterData.h_view().setBubbleAvHeId(this->_numDOFs + 1);
	this->_clusterData.h_view().setBubbleAvVId(this->_numDOFs + 2);
	this->_numDOFs += 3;
}

template <typename TSpeciesEnum>
void
PSIReactionNetwork<TSpeciesEnum>::updateExtraClusterData(
	const std::vector<double>& gridTemps, const std::vector<double>& gridDepths)
{
	if (!this->_enableTrapMutation) {
		return;
	}

	// Check which temperature index to use
	IdType tempId = 0;
	for (tempId = 0; tempId < gridDepths.size(); tempId++) {
		if (gridDepths[tempId] > 0.01)
			break;
	}

	_tmHandler->updateData(gridTemps[tempId]);

	auto& tmData = this->_clusterData.h_view().extraData.trapMutationData;

	using Kokkos::HostSpace;
	using Kokkos::MemoryUnmanaged;

	auto desorpInit = _tmHandler->getDesorptionInitializer();
	auto subpaving = this->_subpaving;
	IndexType desorpId = this->invalidIndex();
	Kokkos::parallel_reduce(
		1,
		KOKKOS_LAMBDA(std::size_t, IndexType & running) {
			Composition comp{};
			comp[Species::He] = desorpInit.size;
			running = static_cast<IndexType>(subpaving.findTileId(comp));
		},
		desorpId);
	auto desorp = create_mirror_view(tmData.desorption);
	desorp() = detail::Desorption{desorpInit, desorpId};
	deep_copy(tmData.desorption, desorp);

	auto depths = Kokkos::View<const double[7], HostSpace, MemoryUnmanaged>(
		_tmHandler->getDepths().data());
	deep_copy(tmData.tmDepths, depths);

	auto vSizes = Kokkos::View<const AmountType[7], HostSpace, MemoryUnmanaged>(
		_tmHandler->getVacancySizes().data());
	deep_copy(tmData.tmVSizes, vSizes);

	this->invalidateDataMirror();
}

template <typename TSpeciesEnum>
void
PSIReactionNetwork<TSpeciesEnum>::selectTrapMutationReactions(
	double depth, double spacing)
{
	auto& tmData = this->_clusterData.h_view().extraData.trapMutationData;
	auto depths = create_mirror_view(tmData.tmDepths);
	deep_copy(depths, tmData.tmDepths);
	auto enable = create_mirror_view(tmData.tmEnabled);
	for (std::size_t l = 0; l < depths.size(); ++l) {
		enable[l] = false;
		if (depths[l] == 0.0) {
			continue;
		}
		if (depths[l] < depth + 0.01 && depths[l] > depth - spacing - 0.01) {
			enable[l] = true;
		}
	}
	deep_copy(tmData.tmEnabled, enable);

	// NOTE:
	// Not calling invalidateDataMirror() here because this change should
	// only matter to reactions on-device
}

template <typename TSpeciesEnum>
void
PSIReactionNetwork<TSpeciesEnum>::computeFluxesPreProcess(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex,
	double surfaceDepth, double spacing)
{
	if (this->_enableTrapMutation) {
		updateDesorptionLeftSideRate(concentrations, gridIndex);
		selectTrapMutationReactions(surfaceDepth, spacing);
	}

	if (this->_enableBursting) {
		this->_clusterData.h_view().setDepth(surfaceDepth);
	}

	if (this->_enableLargeBubble) {
		auto clusterDataMirror = this->getClusterDataMirror();

		// Get the concentrations on the host
		auto dConcs = Kokkos::subview(concentrations,
			std::make_pair(clusterDataMirror.bubbleId(),
				clusterDataMirror.bubbleAvVId() + 1));
		auto hConcs = create_mirror_view(dConcs);
		deep_copy(hConcs, dConcs);

		auto Cb = hConcs(0);
		auto avHe = hConcs(1) / Cb;
		auto avV = hConcs(2) / Cb;
		if (Cb == 0.0) {
			avHe = 0.0;
			avV = 0.0;
		}

		this->_clusterData.h_view().setBubbleAvHe(avHe);
		this->_clusterData.h_view().setBubbleAvV(avV);
		this->_clusterData.h_view().setBubbleAvRadius(util::max(0.0,
			computeBubbleRadius(avV, clusterDataMirror.latticeParameter())));
	}
}

template <typename TSpeciesEnum>
void
PSIReactionNetwork<TSpeciesEnum>::computePartialsPreProcess(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex, double surfaceDepth, double spacing)
{
	if (this->_enableTrapMutation) {
		updateDesorptionLeftSideRate(concentrations, gridIndex);
		selectTrapMutationReactions(surfaceDepth, spacing);
	}

	if (this->_enableBursting) {
		this->_clusterData.h_view().setDepth(surfaceDepth);
	}

	if (this->_enableLargeBubble) {
		auto clusterDataMirror = this->getClusterDataMirror();

		// Get the concentrations on the host
		auto dConcs = Kokkos::subview(concentrations,
			std::make_pair(clusterDataMirror.bubbleId(),
				clusterDataMirror.bubbleAvVId() + 1));
		auto hConcs = create_mirror_view(dConcs);
		deep_copy(hConcs, dConcs);

		auto Cb = hConcs(0);
		auto avHe = hConcs(1) / Cb;
		auto avV = hConcs(2) / Cb;
		if (Cb == 0.0) {
			avHe = 0.0;
			avV = 0.0;
		}

		this->_clusterData.h_view().setBubbleAvHe(avHe);
		this->_clusterData.h_view().setBubbleAvV(avV);
		this->_clusterData.h_view().setBubbleAvRadius(util::max(0.0,
			computeBubbleRadius(avV, clusterDataMirror.latticeParameter())));
	}
}

template <typename TSpeciesEnum>
void
PSIReactionNetwork<TSpeciesEnum>::updateReactionRates(double time)
{
	Superclass::updateReactionRates(time);

	using TrapMutationReactionType =
		typename Superclass::Traits::TrapMutationReactionType;
	// TODO: is this just the local largest rate? Is it correct?
	auto largestRate = this->getLargestRate();
	auto tmReactions =
		this->_reactions.template getView<TrapMutationReactionType>();
	Kokkos::parallel_for(
		"PSIReactionNetwork::updateReactionRates", tmReactions.size(),
		KOKKOS_LAMBDA(
			IndexType i) { tmReactions[i].updateLargestRates(largestRate); });
}

template <typename TSpeciesEnum>
void
PSIReactionNetwork<TSpeciesEnum>::updateTrapMutationDisappearingRate(
	double totalTrappedHeliumConc)
{
	// Set the rate to have an exponential decrease
	if (this->_enableAttenuation) {
		auto& tmData = this->_clusterData.h_view().extraData.trapMutationData;
		auto mirror = create_mirror_view(tmData.currentDisappearingRate);
		mirror() = exp(-4.0 * totalTrappedHeliumConc);
		deep_copy(tmData.currentDisappearingRate, mirror);

		// NOTE:
		// Not calling invalidateDataMirror() here because this change should
		// only matter to reactions on-device
	}
}

template <typename TSpeciesEnum>
void
PSIReactionNetwork<TSpeciesEnum>::updateDesorptionLeftSideRate(
	ConcentrationsView concentrations, IndexType gridIndex)
{
	// TODO: Desorption is constant. So make it available on both host and
	// device. Either DualView or just direct value type that gets copied
	auto& tmData = this->_clusterData.h_view().extraData.trapMutationData;
	auto desorp = create_mirror_view(tmData.desorption);
	deep_copy(desorp, tmData.desorption);
	auto lsRate = create_mirror_view(tmData.currentDesorpLeftSideRate);
	lsRate() = this->getLeftSideRate(concentrations, desorp().id, gridIndex);
	deep_copy(tmData.currentDesorpLeftSideRate, lsRate);

	// NOTE:
	// Not calling invalidateDataMirror() here because this change should
	// only matter to reactions on-device
}

template <typename TSpeciesEnum>
double
PSIReactionNetwork<TSpeciesEnum>::checkLatticeParameter(double latticeParameter)
{
	if (latticeParameter <= 0.0) {
		return tungstenLatticeConstant;
	}
	return latticeParameter;
}

template <typename TSpeciesEnum>
double
PSIReactionNetwork<TSpeciesEnum>::checkImpurityRadius(double impurityRadius)
{
	if (impurityRadius <= 0.0) {
		return heliumRadius;
	}
	return impurityRadius;
}

template <typename TSpeciesEnum>
typename PSIReactionNetwork<TSpeciesEnum>::IndexType
PSIReactionNetwork<TSpeciesEnum>::checkLargestClusterId()
{
	// Copy the cluster data for the parallel loop
	auto clData = this->_clusterData.d_view;
	using Reducer = Kokkos::MaxLoc<PSIReactionNetwork<TSpeciesEnum>::AmountType,
		PSIReactionNetwork<TSpeciesEnum>::IndexType>;
	typename Reducer::value_type maxLoc;
	Kokkos::parallel_reduce(
		"PSIReactionNetwork::checkLargestClusterId", this->_numClusters,
		KOKKOS_LAMBDA(IndexType i, typename Reducer::value_type & update) {
			const auto& clReg = clData().getCluster(i).getRegion();
			Composition hi = clReg.getUpperLimitPoint();
			auto size = hi[Species::He] + hi[Species::V];
			if constexpr (psi::hasDeuterium<Species>) {
				size += hi[Species::D];
			}
			if constexpr (psi::hasTritium<Species>) {
				size += hi[Species::T];
			}
			if (size > update.val) {
				update.val = size;
				update.loc = i;
			}
		},
		Reducer(maxLoc));

	return maxLoc.loc;
}

namespace detail
{
template <typename TSpeciesEnum>
PSIReactionGenerator<TSpeciesEnum>::PSIReactionGenerator(
	const PSIReactionNetwork<TSpeciesEnum>& network) :
	Superclass(network),
	largestClusterId(network.largestClusterId)
{
	bool enableTrapMutation = network.getEnableTrapMutation();
	if (!enableTrapMutation) {
		return;
	}

	auto tmHandler =
		psi::getTrapMutationHandler(enableTrapMutation, network.getMaterial());
	auto vSizes = tmHandler->getAllVSizes();
	for (std::size_t n = 0; n < vSizes.size(); ++n) {
		if (vSizes[n].empty()) {
			continue;
		}
		auto nv = vSizes[n].size();
		_tmVSizes[n] = Kokkos::View<AmountType*>(
			"TM Vacancy Sizes - He" + std::to_string(n + 1), nv);
		auto vSizesMirror = Kokkos::View<AmountType*, Kokkos::HostSpace,
			Kokkos::MemoryUnmanaged>(vSizes[n].data(), nv);
		deep_copy(_tmVSizes[n], vSizesMirror);
	}
}

template <typename TSpeciesEnum>
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
PSIReactionGenerator<TSpeciesEnum>::operator()(
	IndexType i, IndexType j, TTag tag) const
{
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;

	constexpr auto species = NetworkType::getSpeciesRange();
	constexpr auto speciesNoI = NetworkType::getSpeciesRangeNoI();

	auto numClusters = this->getNumberOfClusters();

	// Check the diffusion factors
	auto diffusionFactor = this->_clusterData.diffusionFactor;
	if (i == j) {
		addBurstings(i, tag);
		if (diffusionFactor(i) > 0.0)
			addSinks(i, tag);
	}

	if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
		return;
	}
	if (this->_clusterData.enableLargeBubble())
		addLargeBubbleReactions(i, j, tag);

	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	auto& subpaving = this->getSubpaving();

	// Special case for I + I
	if (lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::I)) {
		// Compute the composition of the new cluster
		auto minSize = lo1[Species::I] + lo2[Species::I];
		auto maxSize = hi1[Species::I] + hi2[Species::I] - 2;
		auto prodBounds = Kokkos::make_pair<int, int>(0, 0);
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			// Check the bounds of the previous product
			if (k >= prodBounds.first && k <= prodBounds.second) {
				// Do nothing because we already added this reaction
				continue;
			}
			Composition comp = Composition::zero();
			comp[Species::I] = k;
			auto iProdId = subpaving.findTileId(comp);
			if (iProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, iProdId});
				if (lo1[Species::I] == 1 || lo2[Species::I] == 1) {
					this->addDissociationReaction(tag, {iProdId, i, j});
				}
				// Update the bounds
				const auto& prodReg = this->getCluster(iProdId).getRegion();
				Composition loP = prodReg.getOrigin();
				Composition hiP = prodReg.getUpperLimitPoint();
				prodBounds.first = loP[Species::I];
				prodBounds.second = hiP[Species::I] - 1;
			}
		}
		return;
	}

	// Special case for I + V
	if ((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::I))) {
		// Find out which one is which
		auto vReg = lo1.isOnAxis(Species::V) ? cl1Reg : cl2Reg;
		auto iReg = lo1.isOnAxis(Species::I) ? cl1Reg : cl2Reg;
		// Compute the largest possible product and the smallest one
		int largestProd =
			(int)vReg[Species::V].end() - 1 - (int)iReg[Species::I].begin();
		int smallestProd =
			(int)vReg[Species::V].begin() - (int)iReg[Species::I].end() + 1;
		auto prodBounds = Kokkos::make_pair<int, int>(0, 0);
		// Loop on the products
		for (int prodSize = smallestProd; prodSize <= largestProd; prodSize++) {
			// Check the bounds of the previous product
			if (prodSize >= prodBounds.first and
				prodSize <= prodBounds.second and prodSize != 0) {
				// Do nothing because we already added this reaction
				continue;
			}
			// 3 cases
			if (prodSize > 0) {
				// Looking for V cluster
				Composition comp = Composition::zero();
				comp[Species::V] = prodSize;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex()) {
					this->addProductionReaction(tag, {i, j, vProdId});
					// No dissociation

					// Update the bounds
					const auto& prodReg = this->getCluster(vProdId).getRegion();
					Composition loP = prodReg.getOrigin();
					Composition hiP = prodReg.getUpperLimitPoint();
					prodBounds.first = loP[Species::V];
					prodBounds.second = hiP[Species::V] - 1;
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

					// Update the bounds
					const auto& prodReg = this->getCluster(iProdId).getRegion();
					Composition loP = prodReg.getOrigin();
					Composition hiP = prodReg.getUpperLimitPoint();
					prodBounds.first = 1 - (int)hiP[Species::I];
					prodBounds.second = -(int)loP[Species::I];
				}
			}
			else {
				// No product
				this->addProductionReaction(tag, {i, j});
				// Update the bounds
				prodBounds.first = 0;
				prodBounds.second = 0;
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
		// TODO: check l correspond to the same species in bounds and prod
		for (auto l : speciesNoI) {
			if (prodReg[l()].begin() > bounds[l].second) {
				isGood = false;
				break;
			}
			if (prodReg[l()].end() - 1 < bounds[l].first) {
				isGood = false;
				break;
			}
		}

		if (isGood) {
			// Increase nProd
			nProd++;
			this->addProductionReaction(tag, {i, j, k});
			// TODO: will have to add some rules, i or j should be a simplex
			// cluster of max size 1
			if (!cl1Reg.isSimplex() || !cl2Reg.isSimplex() ||
				!prodReg.isSimplex()) {
				continue;
			}
			// Loop on the species
			bool isOnAxis1 = false, isOnAxis2 = false;
			for (auto l : species) {
				if (lo1.isOnAxis(l()) && lo1[l()] == 1)
					isOnAxis1 = true;
				if (lo2.isOnAxis(l()) && lo2[l()] == 1)
					isOnAxis2 = true;
			}
			if (isOnAxis1 || isOnAxis2) {
				if (lo1.isOnAxis(Species::I) || lo2.isOnAxis(Species::I)) {
					continue;
				}

				this->addDissociationReaction(tag, {k, i, j});
			}
		}
	}

	// Modified Trap-Mutation
	if (this->_clusterData.enableTrapMutation()) {
		auto heAmt = lo1[Species::He];
		if (cl1Reg.isSimplex() && cl2Reg.isSimplex() && 1 <= heAmt &&
			heAmt <= 7) {
			Composition comp1 = Composition::zero();
			comp1[Species::He] = heAmt;
			auto comp2 = comp1;
			auto& vSizes = _tmVSizes[heAmt - 1];
			for (std::size_t n = 0; n < vSizes.extent(0); ++n) {
				auto& vSize = vSizes[n];
				comp2[Species::V] = vSize;
				if (lo1 == comp1 && lo2 == comp2) {
					Composition compI = Composition::zero();
					compI[Species::I] = vSize;
					auto iClusterId = subpaving.findTileId(compI);
					this->addTrapMutationReaction(tag, {i, j, iClusterId});
				}
			}
		}
	}

	// Special case for trap-mutation
	if (nProd == 0) {
		// Look for larger clusters only if one of the reactant is pure He
		if (!(cl1Reg.isSimplex() && lo1.isOnAxis(Species::He)) &&
			!(cl2Reg.isSimplex() && lo2.isOnAxis(Species::He))) {
			return;
		}

		// Check that both reactants contain He
		if (cl1Reg[Species::He].begin() < 1 ||
			cl2Reg[Species::He].begin() < 1) {
			return;
		}

		// Loop on possible I sizes
		// TODO: get the correct value for maxISize
		AmountType maxISize = 6;
		for (AmountType n = 1; n <= maxISize; ++n) {
			// Find the corresponding cluster
			Composition comp = Composition::zero();
			comp[Species::I] = n;
			auto iClusterId = subpaving.findTileId(comp);

			// Check the I cluster exists
			if (iClusterId == NetworkType::invalidIndex())
				continue;

			bounds[Species::V].first += 1;
			bounds[Species::V].second += 1;

			// Look for potential product
			IndexType nProd = 0;
			for (IndexType k = 0; k < numClusters; ++k) {
				// Get the composition
				const auto& prodReg = this->getCluster(k).getRegion();
				bool isGood = true;
				// Loop on the species
				// TODO: check l correspond to the same species in bounds
				// and prod
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
					this->addProductionReaction(tag, {i, j, k, iClusterId});
					// No dissociation
				}
			}
			// Stop if we found a product
			if (nProd > 0) {
				break;
			}
		}
	}
}

template <typename TSpeciesEnum>
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
PSIReactionGenerator<TSpeciesEnum>::addSinks(IndexType i, TTag tag) const
{
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	const auto& clReg = this->getCluster(i).getRegion();
	Composition lo = clReg.getOrigin();

	// I
	if (lo.isOnAxis(Species::I)) {
		this->addSinkReaction(tag, {i, NetworkType::invalidIndex()});
	}

	// V
	if (clReg.isSimplex() && lo.isOnAxis(Species::V)) {
		if (lo[Species::V] == 1)
			this->addSinkReaction(tag, {i, NetworkType::invalidIndex()});
	}
}

template <typename TSpeciesEnum>
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
PSIReactionGenerator<TSpeciesEnum>::addLargeBubbleReactions(
	IndexType i, IndexType j, TTag tag) const
{
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	auto bubbleId = this->_clusterData.bubbleId();

	if (i == j) {
		const auto& clReg = this->getCluster(i).getRegion();
		Composition lo = clReg.getOrigin();

		// Check reaction with largest bubble
		if (not clReg.isSimplex())
			return;

		// He case
		if (lo.isOnAxis(Species::He)) {
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
		// V case
		else if (lo.isOnAxis(Species::V)) {
			// V_k + B -> B
			this->addProductionReaction(tag, {i, bubbleId, bubbleId});
		}
		//		// I case
		//		else if (lo.isOnAxis(Species::I)) {
		//			// I_k + B -> B
		//			this->addProductionReaction(tag, {i, bubbleId, bubbleId});
		//		}
	}

	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	// Special case for I + I
	if (lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::I))
		return;

	// Special case for I + V
	if ((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::I)))
		return;

	// Find the edge of the phase space
	const auto& largestReg =
		this->getCluster(this->largestClusterId).getRegion();
	Composition hiLargest = largestReg.getUpperLimitPoint();
	auto largestV = hiLargest[Species::V] - 1;
	auto largestHe =
		util::getMaxHePerVLoop(largestV, this->_clusterData.latticeParameter(),
			this->_clusterData.getTemperature());

	// General case
	constexpr auto species = NetworkType::getSpeciesRange();
	constexpr auto numSpeciesNoI = NetworkType::getNumberOfSpeciesNoI();
	using BoundsArray = Kokkos::Array<Kokkos::pair<int, int>, numSpeciesNoI>;
	plsm::EnumIndexed<BoundsArray, Species> bounds;
	// Loop on the species
	for (auto l : species) {
		auto low = lo1[l] + lo2[l];
		auto high = hi1[l] + hi2[l] - 2;
		// Special case for I
		if (l == Species::I) {
			bounds[Species::V].first -= (int)high;
			bounds[Species::V].second -= (int)low;
		}
		else {
			bounds[l] = {low, high};
		}
	}

	// If the products is larger than the current phasespace
	if (bounds[Species::He].second > (int)largestHe or
		bounds[Species::V].second > (int)largestV) {
		// Special case for He
		if (lo1.isOnAxis(Species::He) or lo2.isOnAxis(Species::He)) {
			// Is it trap mutation?
			if (bounds[Species::He].first >
				util::getMaxHePerVLoop(bounds[Species::V].first,
					this->_clusterData.latticeParameter(),
					this->getCluster(i).getTemperature(0))) {
				AmountType iSize = 1;
				while (bounds[Species::He].first >
					util::getMaxHePerVLoop(bounds[Species::V].first + iSize,
						this->_clusterData.latticeParameter(),
						this->getCluster(i).getTemperature(0))) {
					iSize++;
				}
				// Get the corresponding I cluster
				auto& subpaving = this->getSubpaving();
				Composition comp = Composition::zero();
				comp[Species::I] = iSize;
				auto iClusterId = subpaving.findTileId(comp);
				if (iClusterId != NetworkType::invalidIndex()) {
					this->addProductionReaction(
						tag, {i, j, bubbleId, iClusterId});
				}
				else {
					this->addProductionReaction(tag, {i, j, bubbleId});
				}
			}
			else {
				// Standard absorption
				this->addProductionReaction(tag, {i, j, bubbleId});
			}
		}
		// Add the reaction
		else {
			this->addProductionReaction(tag, {i, j, bubbleId});
		}
	}

	//	// Special case for I
	//	if (lo1.isOnAxis(Species::I) or lo2.isOnAxis(Species::I)) {
	//		// Find out which one is which
	//		auto hevLo = lo1.isOnAxis(Species::I) ? lo2 : lo1;
	//		auto iLo = lo1.isOnAxis(Species::I) ? lo1 : lo2;
	//
	//		if (iLo[Species::I] + hevLo[Species::V] > largestV) {
	//			auto hevId = lo1.isOnAxis(Species::I) ? j : i;
	//			auto iId = lo1.isOnAxis(Species::I) ? i : j;
	//
	//			// Add the reaction
	//			this->addProductionReaction(tag, {iId, bubbleId, hevId});
	//		}
	//	}
}

template <typename TSpeciesEnum>
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
PSIReactionGenerator<TSpeciesEnum>::addBurstings(IndexType i, TTag tag) const
{
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Large bubble, add it once
	if (i == 0) {
		auto bubbleId = this->_clusterData.bubbleId();
		this->addBurstingReaction(tag, {bubbleId, bubbleId});
	}

	const auto& clReg = this->getCluster(i).getRegion();
	Composition lo = clReg.getOrigin();
	Composition hi = clReg.getUpperLimitPoint();

	// Need helium
	if (hi[Species::He] == 1)
		return;

	// Loop on V
	auto previousIndex = NetworkType::invalidIndex();
	for (auto nV = lo[Species::V]; nV < hi[Species::V]; nV++) {
		// Pure helium case
		if (nV == 0) {
			return;
		}
		// Bubble case
		else {
			auto& subpaving = this->getSubpaving();
			// Look for the V cluster of the same size
			Composition comp = Composition::zero();
			comp[Species::V] = nV;
			auto vClusterId = subpaving.findTileId(comp);
			if (vClusterId != NetworkType::invalidIndex() and
				vClusterId != previousIndex) {
				this->addBurstingReaction(tag, {i, vClusterId});
				previousIndex = vClusterId;
			}
		}
	}
}

template <typename TSpeciesEnum>
inline ReactionCollection<
	typename PSIReactionGenerator<TSpeciesEnum>::NetworkType>
PSIReactionGenerator<TSpeciesEnum>::getReactionCollection() const
{
	ReactionCollection<NetworkType> ret(this->_clusterData.gridSize,
		this->_clusterData.numClusters, this->_enableReadRates,
		this->getProductionReactions(), this->getDissociationReactions(),
		this->getSinkReactions(), this->getTrapMutationReactions(),
		this->getBurstingReactions());
	return ret;
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
