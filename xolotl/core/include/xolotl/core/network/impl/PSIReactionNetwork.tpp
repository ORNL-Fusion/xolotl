#pragma once

#include <xolotl/core/network/detail/impl/SinkReactionGenerator.tpp>
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
void
PSIReactionNetwork<TSpeciesEnum>::updateBurstingConcs(
	double* gridPointSolution, double factor, std::vector<double>& nBurst)
{
	using detail::toIndex;

	// Loop on every cluster
	for (unsigned int i = 0; i < this->getNumClusters(); i++) {
		const auto& clReg = this->getCluster(i, plsm::onHost).getRegion();
		// Non-grouped clusters
		if (clReg.isSimplex()) {
			// Get the composition
			Composition comp = clReg.getOrigin();
			// Pure He, D, or T case
			if (comp.isOnAxis(Species::He)) {
				// Compute the number of atoms released
				nBurst[toIndex(Species::He)] +=
					gridPointSolution[i] * (double)comp[Species::He] * factor;
				// Reset concentration
				gridPointSolution[i] = 0.0;
				continue;
			}
			if constexpr (psi::hasDeuterium<Species>) {
				if (comp.isOnAxis(Species::D)) {
					// Compute the number of atoms released
					nBurst[toIndex(Species::D)] += gridPointSolution[i] *
						(double)comp[Species::D] * factor;
					// Reset concentration
					gridPointSolution[i] = 0.0;
					continue;
				}
			}
			if constexpr (psi::hasTritium<Species>) {
				if (comp.isOnAxis(Species::T)) {
					// Compute the number of atoms released
					nBurst[toIndex(Species::T)] += gridPointSolution[i] *
						(double)comp[Species::T] * factor;
					// Reset concentration
					gridPointSolution[i] = 0.0;
					continue;
				}
			}
			// Mixed cluster case
			if (!comp.isOnAxis(Species::V) && !comp.isOnAxis(Species::I)) {
				// Compute the number of atoms released
				nBurst[toIndex(Species::He)] +=
					gridPointSolution[i] * (double)comp[Species::He] * factor;
				if constexpr (psi::hasDeuterium<Species>) {
					nBurst[toIndex(Species::D)] += gridPointSolution[i] *
						(double)comp[Species::D] * factor;
				}
				if constexpr (psi::hasTritium<Species>) {
					nBurst[toIndex(Species::T)] += gridPointSolution[i] *
						(double)comp[Species::T] * factor;
				}
				// Transfer concentration to V of the same size
				Composition vComp = Composition::zero();
				vComp[Species::V] = comp[Species::V];
				auto vCluster = this->findCluster(vComp, plsm::onHost);
				gridPointSolution[vCluster.getId()] += gridPointSolution[i];
				gridPointSolution[i] = 0.0;

				continue;
			}
		}
		// Grouped clusters
		else {
			// Compute the number of atoms released
			double concFactor = clReg.volume() / clReg[Species::He].length();
			for (auto j : makeIntervalRange(clReg[Species::He])) {
				nBurst[toIndex(Species::He)] +=
					gridPointSolution[i] * (double)j * concFactor * factor;
			}
			if constexpr (psi::hasDeuterium<Species>) {
				concFactor = clReg.volume() / clReg[Species::D].length();
				for (auto j : makeIntervalRange(clReg[Species::D])) {
					nBurst[toIndex(Species::D)] +=
						gridPointSolution[i] * (double)j * concFactor * factor;
				}
			}
			if constexpr (psi::hasTritium<Species>) {
				concFactor = clReg.volume() / clReg[Species::T].length();
				for (auto j : makeIntervalRange(clReg[Species::T])) {
					nBurst[toIndex(Species::T)] +=
						gridPointSolution[i] * (double)j * concFactor * factor;
				}
			}

			// Get the factor
			concFactor = clReg.volume() / clReg[Species::V].length();
			// Loop on the Vs
			for (auto j : makeIntervalRange(clReg[Species::V])) {
				// Transfer concentration to V of the same size
				Composition vComp = Composition::zero();
				vComp[Species::V] = j;
				auto vCluster = this->findCluster(vComp, plsm::onHost);
				// TODO: refine formula with V moment
				gridPointSolution[vCluster.getId()] +=
					gridPointSolution[i] * concFactor;
			}

			// Reset the concentration and moments
			gridPointSolution[i] = 0.0;
			auto momentIds = this->getCluster(i, plsm::onHost).getMomentIds();
			for (std::size_t j = 0; j < momentIds.extent(0); j++) {
				gridPointSolution[momentIds(j)] = 0.0;
			}
		}
	}
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
	auto clData = typename PSIReactionNetwork<TSpeciesEnum>::ClusterDataRef(
		this->_clusterData);
	using Reducer = Kokkos::MaxLoc<PSIReactionNetwork<TSpeciesEnum>::AmountType,
		PSIReactionNetwork<TSpeciesEnum>::IndexType>;
	typename Reducer::value_type maxLoc;
	Kokkos::parallel_reduce(
		this->_numClusters,
		KOKKOS_LAMBDA(IndexType i, typename Reducer::value_type & update) {
			const auto& clReg = clData.getCluster(i).getRegion();
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

	if (i == j) {
		addSinks(i, tag);
	}

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
	if (lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::I)) {
		// Compute the composition of the new cluster
		auto minSize = lo1[Species::I] + lo2[Species::I];
		auto maxSize = hi1[Species::I] + hi2[Species::I] - 2;
		// Find the corresponding clusters
		for (auto k = minSize; k <= maxSize; k++) {
			Composition comp = Composition::zero();
			comp[Species::I] = k;
			auto iProdId = subpaving.findTileId(comp, plsm::onDevice);
			if (iProdId != subpaving.invalidIndex() &&
				iProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, iProdId});
				if (lo1[Species::I] == 1 || lo2[Species::I] == 1) {
					this->addDissociationReaction(tag, {iProdId, i, j});
				}
				previousIndex = iProdId;
			}
		}
		return;
	}

	// Special case for I + V
	if ((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::I))) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		// Int can be grouped
		auto iReg = lo1.isOnAxis(Species::I) ? cl1Reg : cl2Reg;
		for (auto k : makeIntervalRange(iReg[Species::I])) {
			// Compute the product size
			int prodSize = vSize - k;
			// 3 cases
			if (prodSize > 0) {
				// Looking for V cluster
				Composition comp = Composition::zero();
				comp[Species::V] = prodSize;
				auto vProdId = subpaving.findTileId(comp, plsm::onDevice);
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
				auto iProdId = subpaving.findTileId(comp, plsm::onDevice);
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
			auto iClusterId = subpaving.findTileId(comp, plsm::onDevice);

			bounds[Species::V].first += 1;
			bounds[Species::V].second += 1;

			// Look for potential product
			IndexType nProd = 0;
			for (IndexType k = 0; k < numClusters; ++k) {
				// Get the composition
				const auto& prodReg = this->getCluster(k).getRegion();
				bool isGood = true;
				// Loop on the species
				// TODO: check l correspond to the same species in bounds and
				// prod
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
inline ReactionCollection<
	typename PSIReactionGenerator<TSpeciesEnum>::NetworkType>
PSIReactionGenerator<TSpeciesEnum>::getReactionCollection() const
{
	ReactionCollection<NetworkType> ret(this->_clusterData.gridSize,
		this->getProductionReactions(), this->getDissociationReactions(),
		this->getSinkReactions());
	return ret;
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
