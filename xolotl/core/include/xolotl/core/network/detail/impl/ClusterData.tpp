#pragma once

#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename MemSpace>
struct MemSpaceLabelHelper
{
	static constexpr bool onDevice =
		(std::is_same_v<MemSpace, plsm::DeviceMemSpace> &&
			!std::is_same_v<plsm::HostMemSpace, plsm::DeviceMemSpace>);
	static constexpr char label = onDevice ? 'D' : 'H';
};

template <typename MemSpace>
constexpr char memSpaceLabel = MemSpaceLabelHelper<MemSpace>::label;

inline std::string
labelString(const char labelChar)
{
	return std::string(" (") + labelChar + std::string(")");
}

template <typename MemSpace>
inline std::string
labelStr()
{
	return labelString(memSpaceLabel<MemSpace>);
}

template <typename MemSpace>
ClusterDataCommon<MemSpace>::ClusterDataCommon(
	IndexType numClusters_, IndexType gridSize_) :
	numClusters(numClusters_),
	gridSize(gridSize_),
	_floatVals("Floating Point Values" + labelStr<MemSpace>()),
	_intVals("Integer Point Values" + labelStr<MemSpace>()),
	_boolVals("Boolean Values" + labelStr<MemSpace>()),
	temperature("Temperature" + labelStr<MemSpace>(), gridSize),
	reactionRadius("Reaction Radius" + labelStr<MemSpace>(), numClusters),
	formationEnergy("Formation Energy" + labelStr<MemSpace>(), numClusters),
	migrationEnergy("Migration Energy" + labelStr<MemSpace>(), numClusters),
	diffusionFactor("Diffusion Factor" + labelStr<MemSpace>(), numClusters),
	diffusionCoefficient(
		"Diffusion Coefficient" + labelStr<MemSpace>(), numClusters, gridSize)
{
}

template <typename MemSpace>
template <typename TClusterDataCommon>
inline void
ClusterDataCommon<MemSpace>::deepCopy(const TClusterDataCommon& data)
{
	deep_copy(_floatVals, data._floatVals);
	deep_copy(_intVals, data._intVals);
	deep_copy(_boolVals, data._boolVals);
	deep_copy(temperature, data.temperature);
	deep_copy(reactionRadius, data.reactionRadius);
	deep_copy(formationEnergy, data.formationEnergy);
	deep_copy(migrationEnergy, data.migrationEnergy);
	deep_copy(diffusionFactor, data.diffusionFactor);
	deep_copy(diffusionCoefficient, data.diffusionCoefficient);
}

template <typename MemSpace>
inline std::uint64_t
ClusterDataCommon<MemSpace>::getDeviceMemorySize() const noexcept
{
	std::uint64_t ret = 0;

	ret += sizeof(numClusters);
	ret += sizeof(gridSize);
	ret += _floatVals.required_allocation_size();
	ret += _intVals.required_allocation_size();
	ret += _boolVals.required_allocation_size();

	ret += temperature.required_allocation_size(temperature.size());
	ret += reactionRadius.required_allocation_size(reactionRadius.size());
	ret += formationEnergy.required_allocation_size(formationEnergy.size());
	ret += migrationEnergy.required_allocation_size(migrationEnergy.size());
	ret += diffusionFactor.required_allocation_size(diffusionFactor.size());
	ret += diffusionCoefficient.required_allocation_size(
		diffusionCoefficient.extent(0), diffusionCoefficient.extent(1));

	return ret;
}

template <typename MemSpace>
inline void
ClusterDataCommon<MemSpace>::setGridSize(IndexType gridSize_)
{
	gridSize = gridSize_;
	temperature = View<double*>("Temperature" + labelStr<MemSpace>(), gridSize);
	diffusionCoefficient = View<double**>(
		"Diffusion Coefficient" + labelStr<MemSpace>(), numClusters, gridSize);
}

template <typename TNetwork, typename MemSpace>
ClusterData<TNetwork, MemSpace>::ClusterData(
	const TilesView& tiles_, IndexType numClusters_, IndexType gridSize_) :
	Superclass(numClusters_, gridSize_),
	tiles(tiles_),
	momentIds(Kokkos::ViewAllocateWithoutInitializing(
				  "Moment Ids" + labelStr<MemSpace>()),
		numClusters_)
{
}

template <typename TNetwork, typename MemSpace>
ClusterData<TNetwork, MemSpace>::ClusterData(
	Subpaving& subpaving, IndexType gridSize_) :
	ClusterData(subpaving.getTiles(), subpaving.getNumberOfTiles(), gridSize_)
{
}

template <typename TNetwork, typename MemSpace>
template <typename TClusterData>
inline void
ClusterData<TNetwork, MemSpace>::deepCopy(const TClusterData& data)
{
	Superclass::deepCopy(data);

	// NOTE: Intentionally omitting tiles assuming that was part of construction
	deep_copy(momentIds, data.momentIds);

	extraData.deepCopy(data.extraData);
}

template <typename TNetwork, typename MemSpace>
inline std::uint64_t
ClusterData<TNetwork, MemSpace>::getDeviceMemorySize() const noexcept
{
	std::uint64_t ret = Superclass::getDeviceMemorySize();

	ret += tiles.required_allocation_size(tiles.size());
	ret += momentIds.required_allocation_size(momentIds.size());
	ret += extraData.getDeviceMemorySize();

	return ret;
}

template <typename TNetwork, typename MemSpace>
inline void
ClusterData<TNetwork, MemSpace>::generate(const ClusterGenerator& generator,
	double latticeParameter, double interstitialBias, double impurityRadius)
{
	auto data = *this;
	Kokkos::parallel_for(
		"ClusterData::generate", this->numClusters,
		KOKKOS_LAMBDA(const IndexType i) {
			auto cluster = data.getCluster(i);
			data.formationEnergy(i) = generator.getFormationEnergy(cluster);
			data.migrationEnergy(i) = generator.getMigrationEnergy(cluster);
			data.diffusionFactor(i) =
				generator.getDiffusionFactor(cluster, latticeParameter);
			data.reactionRadius(i) = generator.getReactionRadius(
				cluster, latticeParameter, interstitialBias, impurityRadius);
		});

	Kokkos::fence();
}

template <typename TNetwork, typename MemSpace>
inline void
ClusterData<TNetwork, MemSpace>::updateDiffusionCoefficients()
{
	using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
	auto data = *this;
	auto updater = ClusterUpdater{};
	Kokkos::parallel_for(
		"ClusterData::updateDiffusionCoefficients",
		Range2D({0, 0}, {this->numClusters, this->gridSize}),
		KOKKOS_LAMBDA(IndexType i, IndexType j) {
			if (!util::equal(data.diffusionFactor(i), 0.0)) {
				updater.updateDiffusionCoefficient(data, i, j);
			}
		});

	Kokkos::fence();
}

template <typename TNetwork, typename MemSpace>
inline typename ClusterData<TNetwork, MemSpace>::IndexType
ClusterData<TNetwork, MemSpace>::defineMomentIds()
{
	constexpr auto speciesRange = TNetwork::getSpeciesRangeForGrouping();

	auto numClusters = this->numClusters;
	auto counts = Kokkos::View<IndexType*>("Moment Id Counts", numClusters);

	auto data = *this;
	IndexType nMomentIds = 0;
	Kokkos::parallel_reduce(
		"ClusterData::defineMomentIds::count", numClusters,
		KOKKOS_LAMBDA(const IndexType i, IndexType& running) {
			const auto& reg = data.getCluster(i).getRegion();
			IndexType count = 0;
			for (auto k : speciesRange) {
				if (reg[k].length() != 1) {
					++count;
				}
			}
			running += count;
			counts(i) = count;
		},
		nMomentIds);

	Kokkos::parallel_scan(
		"ClusterData::defineMomentIds::scan", numClusters,
		KOKKOS_LAMBDA(IndexType i, IndexType & update, const bool finalPass) {
			const auto temp = counts(i);
			if (finalPass) {
				counts(i) = update;
			}
			update += temp;
		});

	Kokkos::parallel_for(
		"ClusterData::defineMomentIds::assignMomentIds", numClusters,
		KOKKOS_LAMBDA(const IndexType i) {
			for (IndexType j = 0; j < data.momentIds.extent(1); ++j) {
				data.momentIds(i, j) = TNetwork::invalidIndex();
			}

			const auto& reg = data.getCluster(i).getRegion();
			IndexType current = counts(i);
			for (auto k : speciesRange) {
				if (reg[k].length() == 1) {
					if (data.momentIds(i, TNetwork::mapToMomentId(k)) ==
						numClusters + current - 1)
						continue;
					data.momentIds(i, TNetwork::mapToMomentId(k)) =
						TNetwork::invalidIndex();
				}
				else {
					data.momentIds(i, TNetwork::mapToMomentId(k)) =
						numClusters + current;
					++current;
				}
			}
		});

	Kokkos::fence();
	return nMomentIds;
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
