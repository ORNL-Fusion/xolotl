#pragma once

#include <type_traits>

#include <Kokkos_Core.hpp>

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#define DEVICE_LAMBDA [=] __device__
#else
#define DEVICE_LAMBDA [=]
#endif

#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/detail/ClusterSet.h>
#include <xolotl/core/network/detail/MultiElementCollection.h>
#include <xolotl/core/network/detail/ReactionData.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename TNetwork>
class ReactionCollection
{
public:
	using NetworkType = TNetwork;
	using IndexType = detail::ReactionNetworkIndexType;
	using ReactionTypes = ReactionTypeList<NetworkType>;
	using Types = ReactionNetworkTypes<NetworkType>;
	using ClusterData = typename Types::ClusterData;
	using RateVector = IReactionNetwork::RateVector;

private:
	static constexpr std::size_t numReactionTypes =
		std::tuple_size<ReactionTypes>::value;

	static constexpr std::size_t numSpeciesNoI =
		NetworkType::getNumberOfSpeciesNoI();

public:
	ReactionCollection() = default;

	template <typename... TViews>
	ReactionCollection(IndexType gridSize, TViews... views) :
		_reactions(views...),
		_data(_reactions.getNumberOfElements(), gridSize,
			_reactions.getElementBeginIndices())
	{
		static_assert(sizeof...(TViews) == numReactionTypes,
			"Construction from views requires the number of views to match the "
			"number of reaction types in the ReactionCollection");

		_reactions.forEachType(
			[this](IndexType reactionTypeIndex, IndexType numReactions,
				auto reactionTypeTag) {
				using ReactionType = typename decltype(reactionTypeTag)::Type;
				_data.coeffs[reactionTypeIndex] =
					ReactionType::allocateCoefficientsView(numReactions);
			});
	}

	void
	setGridSize(IndexType gridSize)
	{
		_data.setGridSize(gridSize);
	}

	void
	setConnectivity(const ClusterConnectivity<>& connectivity)
	{
		auto conn = connectivity;
		forEach(
			"ReactionCollection::setConnectivity",
			DEVICE_LAMBDA(
				auto&& reaction) { reaction.defineJacobianEntries(conn); });
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		std::uint64_t ret = _reactions.getDeviceMemorySize();
		ret += _data.getDeviceMemorySize();
		return ret;
	}

	IndexType
	getNumberOfReactions() const noexcept
	{
		return _data.numReactions;
	}

	template <typename TReaction>
	Kokkos::View<TReaction*>
	getView() const
	{
		return _reactions.template getView<TReaction>();
	}

	template <typename TReaction>
	void
	setView(Kokkos::View<TReaction*> view)
	{
		_reactions.setView(view);
		_data.numReactions = _reactions.getNumberOfElements();
	}

	void
	constructAll(Kokkos::View<ClusterData> clusterData,
		Kokkos::View<ClusterSet*> clusterSets)
	{
		auto chain = _reactions.getChain();
		auto reactionData = ReactionDataRef<NetworkType>(_data);
		// TODO: Enable this without getting the chain
		Kokkos::parallel_for(
			"ReactionCollection::constructAll", _data.numReactions,
			DEVICE_LAMBDA(const IndexType i) {
				chain.apply(
					DEVICE_LAMBDA(auto& reaction) {
						using ReactionType =
							std::remove_reference_t<decltype(reaction)>;
						reaction = ReactionType(
							reactionData, clusterData(), i, clusterSets(i));
					},
					i);
			});
	}

	void
	updateAll(Kokkos::View<ClusterData> clusterData)
	{
		auto reactionData = ReactionDataRef<NetworkType>(_data);
		forEach(
			"ReactionCollection::clusterData", DEVICE_LAMBDA(auto&& reaction) {
				reaction.updateData(reactionData, clusterData());
			});
	}

	void
	updateRates(double time = 0.0)
	{
		forEach(
			"ReactionNetwork::updateReactionRates",
			DEVICE_LAMBDA(auto&& reaction) { reaction.updateRates(time); });
		Kokkos::fence();
	}

	double
	getLargestRate() const
	{
		using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
		double largestRate = 0.0;
		auto rates = _data.rates;
		auto nRates = rates.extent(0);
		auto gridSize = rates.extent(1);
		Kokkos::parallel_reduce(
			"ReactionCollection::getLargestRate",
			Range2D({0, 0}, {nRates, gridSize}),
			KOKKOS_LAMBDA(const IndexType i, const IndexType j, double& max) {
				if (rates(i, j) > max) {
					max = rates(i, j);
				}
			},
			Kokkos::Max<double>(largestRate));
		Kokkos::fence();
		return largestRate;
	}

	template <typename F>
	void
	forEach(const F& func)
	{
		_reactions.forEach(func);
	}

	template <typename F>
	void
	forEach(const std::string& label, const F& func)
	{
		_reactions.forEach(label, func);
	}

	template <typename TReaction, typename F>
	void
	forEachOn(const F& func)
	{
		_reactions.template forEachOn<TReaction>(func);
	}

	template <typename TReaction, typename F>
	void
	forEachOn(const std::string& label, const F& func)
	{
		_reactions.template forEachOn<TReaction>(label, func);
	}

	template <typename F, typename T>
	void
	reduce(const F& func, T& out)
	{
		_reactions.reduce(func, out);
	}

	template <typename F, typename T>
	void
	reduce(const std::string& label, const F& func, T& out)
	{
		_reactions.reduce(label, func, out);
	}

	template <typename TReaction, typename F, typename T>
	void
	reduceOn(const F& func, T& out)
	{
		_reactions.template reduceOn<TReaction>(func, out);
	}

	template <typename TReaction, typename F, typename T>
	void
	reduceOn(const std::string& label, const F& func, T& out)
	{
		_reactions.template reduceOn<TReaction>(label, func, out);
	}

private:
	MultiElementCollection<ReactionTypes> _reactions;
	ReactionData<NetworkType> _data;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
