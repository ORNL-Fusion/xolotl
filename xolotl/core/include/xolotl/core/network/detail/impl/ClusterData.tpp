#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename PlsmContext>
struct ContextLabelHelper;

template <>
struct ContextLabelHelper<plsm::OnHost>
{
	static constexpr char label = 'H';
};

template <>
struct ContextLabelHelper<plsm::OnDevice>
{
	static constexpr char label = 'D';
};

template <typename PlsmContext>
constexpr char contextLabel = ContextLabelHelper<PlsmContext>::label;

inline std::string
labelString(const char labelChar)
{
	return std::string(" (") + labelChar + std::string(")");
}

template <typename PlsmContext>
inline std::string
labelStr()
{
	return labelString(contextLabel<PlsmContext>);
}

template <typename PlsmContext, template <typename> typename ViewConvert>
ClusterDataCommon<PlsmContext, ViewConvert>::ClusterDataCommon(
	IndexType numClusters_, IndexType gridSize_) :
	numClusters(numClusters_),
	gridSize(gridSize_),
	_floatVals("Floating Point Values" + labelStr<PlsmContext>()),
	_boolVals("Boolean Values" + labelStr<PlsmContext>()),
	temperature("Temperature" + labelStr<PlsmContext>(), gridSize),
	reactionRadius("Reaction Radius" + labelStr<PlsmContext>(), numClusters),
	formationEnergy("Formation Energy" + labelStr<PlsmContext>(), numClusters),
	migrationEnergy("Migration Energy" + labelStr<PlsmContext>(), numClusters),
	diffusionFactor("Diffusion Factor" + labelStr<PlsmContext>(), numClusters),
	diffusionCoefficient("Diffusion Coefficient" + labelStr<PlsmContext>(),
		numClusters, gridSize)
{
}

template <typename PlsmContext, template <typename> typename ViewConvert>
template <typename TClusterDataCommon>
inline void
ClusterDataCommon<PlsmContext, ViewConvert>::deepCopy(
	const TClusterDataCommon& data)
{
	deep_copy(_floatVals, data._floatVals);
	deep_copy(_boolVals, data._boolVals);
	deep_copy(temperature, data.temperature);
	deep_copy(reactionRadius, data.reactionRadius);
	deep_copy(formationEnergy, data.formationEnergy);
	deep_copy(migrationEnergy, data.migrationEnergy);
	deep_copy(diffusionFactor, data.diffusionFactor);
	deep_copy(diffusionCoefficient, data.diffusionCoefficient);
}

template <typename PlsmContext, template <typename> typename ViewConvert>
inline std::uint64_t
ClusterDataCommon<PlsmContext, ViewConvert>::getDeviceMemorySize()
	const noexcept
{
	std::uint64_t ret = 0;

	ret += sizeof(numClusters);
	ret += sizeof(gridSize);
	ret += _floatVals.required_allocation_size();
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

template <typename PlsmContext, template <typename> typename ViewConvert>
inline void
ClusterDataCommon<PlsmContext, ViewConvert>::setGridSize(IndexType gridSize_)
{
	gridSize = gridSize_;
	temperature =
		View<double*>("Temperature" + labelStr<PlsmContext>(), gridSize);
	diffusionCoefficient =
		View<double**>("Diffusion Coefficient" + labelStr<PlsmContext>(),
			numClusters, gridSize);
}

template <typename TNetwork, typename PlsmContext,
	template <typename> typename ViewConvert>
ClusterDataImpl<TNetwork, PlsmContext, ViewConvert>::ClusterDataImpl(
	const TilesView& tiles_, IndexType numClusters_, IndexType gridSize_) :
	Superclass(numClusters_, gridSize_),
	tiles(tiles_),
	momentIds("Moment Ids" + labelStr<PlsmContext>(), numClusters_)
{
}

template <typename TNetwork, typename PlsmContext,
	template <typename> typename ViewConvert>
ClusterDataImpl<TNetwork, PlsmContext, ViewConvert>::ClusterDataImpl(
	Subpaving& subpaving, IndexType gridSize_) :
	ClusterDataImpl(subpaving.getTiles(PlsmContext{}),
		subpaving.getNumberOfTiles(PlsmContext{}), gridSize_)
{
}

template <typename TNetwork, typename PlsmContext,
	template <typename> typename ViewConvert>
template <typename TClusterData>
inline void
ClusterDataImpl<TNetwork, PlsmContext, ViewConvert>::deepCopy(
	const TClusterData& data)
{
	Superclass::deepCopy(data);

	// NOTE: Intentionally omitting tiles assuming that was part of construction
	deep_copy(momentIds, data.momentIds);

	extraData.deepCopy(data.extraData);
}

template <typename TNetwork, typename PlsmContext,
	template <typename> typename ViewConvert>
inline std::uint64_t
ClusterDataImpl<TNetwork, PlsmContext, ViewConvert>::getDeviceMemorySize()
	const noexcept
{
	std::uint64_t ret = Superclass::getDeviceMemorySize();

	ret += tiles.required_allocation_size(tiles.size());
	ret += momentIds.required_allocation_size(momentIds.size());
	ret += extraData.getDeviceMemorySize();

	return ret;
}

template <typename TNetwork, typename PlsmContext,
	template <typename> typename ViewConvert>
inline void
ClusterDataImpl<TNetwork, PlsmContext, ViewConvert>::generate(
	const ClusterGenerator& generator, double latticeParameter,
	double interstitialBias, double impurityRadius)
{
	auto data = ClusterDataRef<TNetwork, PlsmContext>(*this);
	Kokkos::parallel_for(
		this->numClusters, KOKKOS_LAMBDA(const IndexType i) {
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
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
