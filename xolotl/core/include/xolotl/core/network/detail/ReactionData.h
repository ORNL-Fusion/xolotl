#pragma once

#include <Kokkos_Core.hpp>

#include <xolotl/core/network/ReactionNetworkTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
using CoefficientsView = Kokkos::View<double*****>;
using CoefficientsViewUnmanaged =
	Kokkos::View<double*****, Kokkos::MemoryUnmanaged>;
using ConstantRateView = Kokkos::View<double****>;
using ConstantRateViewUnmanaged =
	Kokkos::View<double****, Kokkos::MemoryUnmanaged>;

/**
 * @brief Stores all the information needed for a reaction
 * (overlap widths, rates, position in the collection, grouping
 * coefficients, etc.)
 *
 * @tparam TNetwork The network type
 */
template <typename TNetwork>
struct ReactionData
{
	using NetworkType = TNetwork;
	using IndexType = ReactionNetworkIndexType;
	using ReactionTypes = ReactionTypeList<NetworkType>;

	static constexpr std::size_t numReactionTypes =
		std::tuple_size<ReactionTypes>::value;

	static constexpr std::size_t numSpeciesNoI =
		NetworkType::getNumberOfSpeciesNoI();

	ReactionData() = default;

	ReactionData(IndexType nReactions, IndexType gridSize,
		const Kokkos::Array<IndexType, numReactionTypes + 1>& rBeginIds) :
		numReactions(nReactions),
		widths("Reaction Widths", numReactions, numSpeciesNoI),
		rates("Reaction Rates", numReactions, gridSize),
		reactionBeginIndices(rBeginIds)
	{
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		std::uint64_t ret = sizeof(numReactions);
		ret += sizeof(reactionBeginIndices);
		ret +=
			widths.required_allocation_size(widths.extent(0), widths.extent(1));
		ret += rates.required_allocation_size(rates.extent(0), rates.extent(1));

		for (std::size_t r = 0; r < numReactionTypes; ++r) {
			ret += coeffs[r].required_allocation_size(coeffs[r].extent(0),
				coeffs[r].extent(1), coeffs[r].extent(2), coeffs[r].extent(3),
				coeffs[r].extent(4));
			ret += constantRates[r].required_allocation_size(
				constantRates[r].extent(0), constantRates[r].extent(1),
				constantRates[r].extent(2));
		}

		return ret;
	}

	void
	setGridSize(IndexType gridSize)
	{
		rates =
			Kokkos::View<double**>("Reaction Rates", numReactions, gridSize);
	}

	IndexType numReactions{};
	Kokkos::View<double**> widths;
	Kokkos::View<double**> rates;
	Kokkos::Array<IndexType, numReactionTypes + 1> reactionBeginIndices;
	Kokkos::Array<CoefficientsView, numReactionTypes> coeffs;
	Kokkos::Array<ConstantRateView, numReactionTypes> constantRates;
};

template <typename TNetwork>
struct ReactionDataRef
{
	using NetworkType = TNetwork;
	using IndexType = ReactionNetworkIndexType;
	using ReactionTypes = ReactionTypeList<NetworkType>;

	static constexpr std::size_t numReactionTypes =
		std::tuple_size<ReactionTypes>::value;

	ReactionDataRef() = default;

	KOKKOS_INLINE_FUNCTION
	ReactionDataRef(const ReactionData<NetworkType>& data) :
		widths(data.widths),
		rates(data.rates),
		reactionBeginIndices(data.reactionBeginIndices)
	{
		for (std::size_t r = 0; r < numReactionTypes; ++r) {
			coeffs[r] = data.coeffs[r];
			constantRates[r] = data.constantRates[r];
		}
	}

	KOKKOS_INLINE_FUNCTION
	auto
	getCoefficients(IndexType reactionId)
	{
		std::size_t r = 0;
		for (; r < numReactionTypes; ++r) {
			if (reactionId < reactionBeginIndices[r + 1]) {
				break;
			}
		}
		assert(r < numReactionTypes);
		return Kokkos::subview(coeffs[r], reactionId - reactionBeginIndices[r],
			Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
	}

	KOKKOS_INLINE_FUNCTION
	auto
	getConstantRates(IndexType reactionId)
	{
		std::size_t r = 0;
		for (; r < numReactionTypes; ++r) {
			if (reactionId < reactionBeginIndices[r + 1]) {
				break;
			}
		}
		assert(r < numReactionTypes);
		return Kokkos::subview(constantRates[r],
			reactionId - reactionBeginIndices[r], Kokkos::ALL, Kokkos::ALL,
			Kokkos::ALL);
	}

	KOKKOS_INLINE_FUNCTION
	auto
	getWidths(IndexType reactionId)
	{
		return Kokkos::subview(widths, reactionId, Kokkos::ALL);
	}

	KOKKOS_INLINE_FUNCTION
	auto
	getRates(IndexType reactionId)
	{
		return Kokkos::subview(rates, reactionId, Kokkos::ALL);
	}

	Kokkos::View<double**, Kokkos::MemoryUnmanaged> widths;
	Kokkos::View<double**, Kokkos::MemoryUnmanaged> rates;
	Kokkos::Array<IndexType, numReactionTypes + 1> reactionBeginIndices;
	Kokkos::Array<CoefficientsViewUnmanaged, numReactionTypes> coeffs;
	Kokkos::Array<ConstantRateViewUnmanaged, numReactionTypes> constantRates;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
