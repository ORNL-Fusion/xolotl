#pragma once

#include <xolotl/core/network/ClusterData.h>

namespace xolotl
{
namespace core
{
namespace network
{
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
		return asDerived()->_data.reactionRadius(_id);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getFormationEnergy()
	{
		return asDerived()->_data.formationEnergy(_id);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getTemperature(IndexType gridIndex)
	{
		return asDerived()->_data.temperature(gridIndex);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getDiffusionCoefficient(IndexType gridIndex)
	{
		return asDerived()->_data.diffusionCoefficient(_id, gridIndex);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getDiffusionFactor()
	{
		return asDerived()->_data.diffusionFactor(_id);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getMigrationEnergy()
	{
		return asDerived()->_data.migrationEnergy(_id);
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
	IndexType _id{detail::InvalidIndex::value};
};

template <typename PlsmContext>
class ClusterCommon : public ClusterBase<ClusterCommon<PlsmContext>>
{
	friend class ClusterBase<ClusterCommon<PlsmContext>>;

public:
	using Superclass = ClusterBase<ClusterCommon<PlsmContext>>;
	using ClusterData = detail::ClusterDataCommon<PlsmContext>;
	using ClusterDataRef = detail::ClusterDataCommonRef<PlsmContext>;
	using IndexType = typename Superclass::IndexType;

	ClusterCommon() = delete;

	KOKKOS_INLINE_FUNCTION
	ClusterCommon(const ClusterData& data, IndexType id) :
		Superclass(id),
		_data{data}
	{
	}

	KOKKOS_INLINE_FUNCTION
	ClusterCommon(const ClusterDataRef& data, IndexType id) :
		Superclass(id),
		_data{data}
	{
	}

private:
	ClusterDataRef _data;
};

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
	using ClusterDataRef = detail::ClusterDataRef<TNetwork, PlsmContext>;
	using IndexType = typename Superclass::IndexType;

	Cluster() = delete;

	KOKKOS_INLINE_FUNCTION
	Cluster(const ClusterData& data, IndexType id) : Superclass(id), _data{data}
	{
	}

	KOKKOS_INLINE_FUNCTION
	Cluster(const ClusterDataRef& data, IndexType id) :
		Superclass(id),
		_data{data}
	{
	}

	KOKKOS_INLINE_FUNCTION
	Region
	getRegion() const
	{
		return _data.tiles(this->getId()).getRegion();
	}

	KOKKOS_INLINE_FUNCTION
	decltype(auto)
	getMomentIds()
	{
		return Kokkos::subview(_data.momentIds, this->getId(), Kokkos::ALL);
	}

private:
	ClusterDataRef _data;
};
} // namespace network
} // namespace core
} // namespace xolotl
