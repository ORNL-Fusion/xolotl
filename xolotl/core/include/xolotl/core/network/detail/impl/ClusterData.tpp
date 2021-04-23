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
	return std::string("_") + labelChar;
}

template <typename PlsmContext>
inline std::string
labelStr()
{
	return labelString(contextLabel<PlsmContext>);
}

template <typename PlsmContext,
	template <typename> typename ViewConvert>
ClusterDataCommon<PlsmContext, ViewConvert>::ClusterDataCommon(
	IndexType numClusters_, IndexType gridSize_) :
	numClusters(numClusters_),
	gridSize(gridSize_),
	atomicVolume("Atomic Volume" + labelStr<PlsmContext>()),
	latticeParameter("Lattice Parameter" + labelStr<PlsmContext>()),
	fissionRate("Fission Rate" + labelStr<PlsmContext>()),
	zeta("Zeta" + labelStr<PlsmContext>()),
	enableStdReaction("Enable Std Reaction" + labelStr<PlsmContext>()),
	enableReSolution("Enable Re-Solution Process" + labelStr<PlsmContext>()),
	enableNucleation("Enable Nucleation Process" + labelStr<PlsmContext>()),
	enableSink("Enable Sink Process" + labelStr<PlsmContext>()),
	enableTrapMutation(
		"Enable Trap Mutation Process" + labelStr<PlsmContext>()),
	temperature("Temperature" + labelStr<PlsmContext>(), gridSize),
	reactionRadius("Reaction Radius" + labelStr<PlsmContext>(), numClusters),
	formationEnergy("Formation Energy" + labelStr<PlsmContext>(), numClusters),
	migrationEnergy("Migration Energy" + labelStr<PlsmContext>(), numClusters),
	diffusionFactor("Diffusion Factor" + labelStr<PlsmContext>(), numClusters),
	diffusionCoefficient("Diffusion Coefficient" + labelStr<PlsmContext>(),
		numClusters, gridSize)
{
}

template <typename PlsmContext,
	template <typename> typename ViewConvert>
inline std::uint64_t
ClusterDataCommon<PlsmContext, ViewConvert>::getDeviceMemorySize()
	const noexcept
{
	std::uint64_t ret = 0;

	ret += sizeof(numClusters);
	ret += sizeof(gridSize);
	ret += atomicVolume.required_allocation_size();
	ret += latticeParameter.required_allocation_size();
	ret += fissionRate.required_allocation_size();
	ret += zeta.required_allocation_size();
	ret += enableStdReaction.required_allocation_size();
	ret += enableReSolution.required_allocation_size();
	ret += enableNucleation.required_allocation_size();
	ret += enableSink.required_allocation_size();
	ret += enableTrapMutation.required_allocation_size();

	ret += temperature.required_allocation_size(temperature.size());
	ret += reactionRadius.required_allocation_size(reactionRadius.size());
	ret += formationEnergy.required_allocation_size(formationEnergy.size());
	ret += migrationEnergy.required_allocation_size(migrationEnergy.size());
	ret += diffusionFactor.required_allocation_size(diffusionFactor.size());
	ret += diffusionCoefficient.required_allocation_size(
		diffusionCoefficient.extent(0), diffusionCoefficient.extent(1));

	return ret;
}

template <typename PlsmContext,
	template <typename> typename ViewConvert>
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

template <typename TN1, typename PC1, template <typename> typename VC1,
	typename TN2, typename PC2, template <typename> typename VC2>
inline void
deepCopy(ClusterDataImpl<TN1, PC1, VC1> to, ClusterDataImpl<TN2, PC2, VC2> from)
{
	deep_copy(to.atomicVolume, from.atomicVolume);
	deep_copy(to.latticeParameter, from.latticeParameter);
	deep_copy(to.fissionRate, from.fissionRate);
	deep_copy(to.zeta, from.zeta);
	deep_copy(to.enableStdReaction, from.enableStdReaction);
	deep_copy(to.enableReSolution, from.enableReSolution);
	deep_copy(to.enableNucleation, from.enableNucleation);
	deep_copy(to.enableSink, from.enableSink);
	deep_copy(to.enableTrapMutation, from.enableTrapMutation);
	deep_copy(to.temperature, from.temperature);
	deep_copy(to.reactionRadius, from.reactionRadius);
	deep_copy(to.formationEnergy, from.formationEnergy);
	deep_copy(to.migrationEnergy, from.migrationEnergy);
	deep_copy(to.diffusionFactor, from.diffusionFactor);
	deep_copy(to.diffusionCoefficient, from.diffusionCoefficient);

	// NOTE: Intentionally omitting tiles assuming that was part of construction
	deep_copy(to.momentIds, from.momentIds);

	deepCopy(to.extraData, from.extraData);
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
