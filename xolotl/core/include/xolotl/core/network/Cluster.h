#pragma once

#include <xolotl/core/network/detail/ClusterData.h>

namespace xolotl
{
namespace core
{
namespace network
{
/**
 * @brief Base class for cluster properties.
 *
 * @tparam TDerived The derived template class
 */
template <typename TDerived>
class ClusterBase
{
public:
	using IndexType = detail::ReactionNetworkIndexType;

	KOKKOS_INLINE_FUNCTION
	ClusterBase(IndexType id) : _id(id)
	{
	}

	KOKKOS_INLINE_FUNCTION
	double
	getReactionRadius()
	{
		return asDerived()->_data->reactionRadius(_id);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getFormationEnergy()
	{
		return asDerived()->_data->formationEnergy(_id);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getTemperature(IndexType gridIndex)
	{
		return asDerived()->_data->temperature(gridIndex);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getDiffusionCoefficient(IndexType gridIndex)
	{
		return asDerived()->_data->diffusionCoefficient(_id, gridIndex);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getDiffusionFactor()
	{
		return asDerived()->_data->diffusionFactor(_id);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getMigrationEnergy()
	{
		return asDerived()->_data->migrationEnergy(_id);
	}

	KOKKOS_INLINE_FUNCTION
	IndexType
	getId() const noexcept
	{
		return _id;
	}

private:
	KOKKOS_INLINE_FUNCTION
	TDerived*
	asDerived() noexcept
	{
		return static_cast<TDerived*>(this);
	}

private:
	IndexType _id{detail::invalidNetworkIndex};
};

/**
 * @brief Class for cluster properties common to every type of clusters.
 *
 * @tparam PlsmContext Host or Device
 */
template <typename PlsmContext>
class ClusterCommon : public ClusterBase<ClusterCommon<PlsmContext>>
{
	friend class ClusterBase<ClusterCommon<PlsmContext>>;

public:
	using Superclass = ClusterBase<ClusterCommon<PlsmContext>>;
	using ClusterData = detail::ClusterDataCommon<PlsmContext>;
	using IndexType = typename Superclass::IndexType;

	ClusterCommon() = delete;

	KOKKOS_INLINE_FUNCTION
	ClusterCommon(const ClusterData* data, IndexType id) :
		Superclass(id),
		_data{data}
	{
	}

private:
	const ClusterData* _data;
};

/**
 * @brief Network specific cluster class.
 *
 * @tparam TNetwork Network type
 * @tparam PlsmContext Host or Device
 */
template <typename TNetwork, typename PlsmContext>
class Cluster : public ClusterBase<Cluster<TNetwork, PlsmContext>>
{
	friend class ClusterBase<Cluster>;

	using Types = detail::ReactionNetworkTypes<TNetwork>;

public:
	using Superclass = ClusterBase<Cluster<TNetwork, PlsmContext>>;
	using Subpaving = typename Types::Subpaving;
	using Region = typename Subpaving::RegionType;
	using ClusterData = detail::ClusterData<TNetwork, PlsmContext>;
	using IndexType = typename Superclass::IndexType;
	using Composition = typename Types::Composition;

	Cluster() = delete;

	KOKKOS_INLINE_FUNCTION
	Cluster(const ClusterData* data, IndexType id) : Superclass(id), _data{data}
	{
	}

	KOKKOS_INLINE_FUNCTION
	const Region&
	getRegion() const
	{
		return this->_data->tiles(this->getId()).getRegion();
	}

	KOKKOS_INLINE_FUNCTION
	Composition
	getOriginComposition() const
	{
		return getRegion().getOrigin();
	}

	KOKKOS_INLINE_FUNCTION
	decltype(auto)
	getMomentIds()
	{
		return Kokkos::subview(_data->momentIds, this->getId(), Kokkos::ALL);
	}

private:
	const ClusterData* _data;
};
} // namespace network
} // namespace core
} // namespace xolotl
