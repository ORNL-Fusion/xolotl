#pragma once

#include <xolotl/core/network/detail/impl/NucleationReactionGenerator.tpp>
#include <xolotl/core/network/detail/impl/ReSolutionReactionGenerator.tpp>
#include <xolotl/core/network/detail/impl/SinkReactionGenerator.tpp>
#include <xolotl/core/network/impl/NEClusterGenerator.tpp>
#include <xolotl/core/network/impl/NEReaction.tpp>
#include <xolotl/core/network/impl/ReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
double
NEReactionNetwork::checkLatticeParameter(double latticeParameter)
{
	if (latticeParameter <= 0.0) {
		return uraniumDioxydeLatticeConstant;
	}
	return latticeParameter;
}

double
NEReactionNetwork::checkImpurityRadius(double impurityRadius)
{
	if (impurityRadius <= 0.0) {
		return xenonRadius;
	}
	return impurityRadius;
}

NEReactionNetwork::IndexType
NEReactionNetwork::checkLargestClusterId()
{
	// Copy the cluster data for the parallel loop
	auto clData = ClusterDataRef(_clusterData);
	using Reducer = Kokkos::MaxLoc<NEReactionNetwork::AmountType,
		NEReactionNetwork::IndexType>;
	Reducer::value_type maxLoc;
	Kokkos::parallel_reduce(
		_numClusters,
		KOKKOS_LAMBDA(IndexType i, Reducer::value_type & update) {
			const Region& clReg = clData.getCluster(i).getRegion();
			Composition hi = clReg.getUpperLimitPoint();
			auto size = hi[Species::Xe] + hi[Species::V];
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
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
NEReactionGenerator::operator()(IndexType i, IndexType j, TTag tag) const
{
	// TODO: re-solution and nucleation need to be added
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;

	constexpr auto species = NetworkType::getSpeciesRange();
	constexpr auto speciesNoI = NetworkType::getSpeciesRangeNoI();
	constexpr auto invalidIndex = NetworkType::invalidIndex();

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

	// Special case for I + V
	if (cl1Reg.isSimplex() && cl2Reg.isSimplex() &&
		((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
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
			auto vProdId = subpaving.findTileId(comp, plsm::onDevice);
			if (vProdId != invalidIndex) {
				this->addProductionReaction(tag, {i, j, vProdId});
				// No dissociation
			}
		}
		else if (prodSize < 0) {
			// Looking for I cluster
			Composition comp = Composition::zero();
			comp[Species::I] = -prodSize;
			auto iProdId = subpaving.findTileId(comp, plsm::onDevice);
			if (iProdId != invalidIndex) {
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

	//	// Special case for Xe_1 + Xe_1V_1
	//	if (cl1Reg.isSimplex() && cl2Reg.isSimplex() &&
	//		((lo1.isOnAxis(Species::Xe) && lo1[Species::Xe] == 1 &&
	//			 lo2[Species::Xe] == 1 && lo2[Species::V] == 1) ||
	//			(lo2.isOnAxis(Species::Xe) && lo2[Species::Xe] == 1 &&
	//				lo1[Species::Xe] == 1 && lo1[Species::V] == 1))) {
	//		return;
	//	}

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
				if (lo1.isOnAxis(Species::Xe) && lo2.isOnAxis(Species::Xe) &&
					lo1[Species::Xe] == 1 && lo2[Species::Xe] == 1) {
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
NEReactionGenerator::addSinks(IndexType i, TTag tag) const
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
	if (clReg.isSimplex() && lo.isOnAxis(Species::V) && lo[Species::V] < 2) {
		this->addSinkReaction(tag, {i, NetworkType::invalidIndex()});
	}
}

inline ReactionCollection<NEReactionGenerator::NetworkType>
NEReactionGenerator::getReactionCollection() const
{
	ReactionCollection<NetworkType> ret(this->_clusterData.gridSize,
		this->getProductionReactions(), this->getDissociationReactions(),
		this->getReSolutionReactions(), this->getNucleationReactions(),
		this->getSinkReactions());
	return ret;
}
} // namespace detail

inline detail::NEReactionGenerator
NEReactionNetwork::getReactionGenerator() const noexcept
{
	return detail::NEReactionGenerator{*this};
}

namespace detail
{
KOKKOS_INLINE_FUNCTION
void
NEClusterUpdater::updateDiffusionCoefficient(
	const ClusterData& data, IndexType clusterId, IndexType gridIndex) const
{
	if (data.migrationEnergy(clusterId) > 0.0) {
		//		using Species = typename NetworkType::Species;
		//		using Composition = typename NetworkType::Composition;

		//		const auto& clReg = data.tiles(clusterId).getRegion();
		//		Composition lo = clReg.getOrigin();
		//		if (clReg.isSimplex() && lo[Species::V] == 1 && lo[Species::Xe]
		//== 1) {
		//			// Intrinsic diffusion
		//			double kernel = -3.04 / (kBoltzmann *
		// data.temperature(gridIndex)); 			double D3 = 7.6e8 *
		// exp(kernel);
		// // nm2/s
		//
		//			// We need the fission rate now
		//			double fissionRate = data.fissionRate(0) * 1.0e27; // #/m3/s
		//
		//			// Athermal diffusion
		//			double D1 = (8e-40 * fissionRate) * 1.0e18; // nm2/s
		//
		//			// Radiation-enhanced diffusion
		//			kernel = -1.2 / (kBoltzmann * data.temperature(gridIndex));
		//			double D2 =
		//				(5.6e-25 * sqrt(fissionRate) * exp(kernel)) * 1.0e18; //
		// nm2/s
		//
		//			data.diffusionCoefficient(clusterId, gridIndex) = D1 + D2 +
		// D3;
		//
		//			return;
		//		}

		data.diffusionCoefficient(clusterId, gridIndex) =
			data.diffusionFactor(clusterId) *
			exp(-data.migrationEnergy(clusterId) /
				(kBoltzmann * data.temperature(gridIndex)));

		return;
	}

	// If the diffusivity is given
	data.diffusionCoefficient(clusterId, gridIndex) =
		data.diffusionFactor(clusterId);
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
