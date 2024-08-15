#pragma once

#include <tuple>

#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/detail/ClusterData.h>
#include <xolotl/core/network/detail/TrapMutationClusterData.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace psi
{
namespace detail
{
template <typename TSpeciesEnum>
struct HasDeuterium : std::false_type
{
};

template <typename TSpeciesEnum>
struct HasTritium : std::false_type
{
};
} // namespace detail

template <typename TSpeciesEnum>
inline constexpr bool hasDeuterium = detail::HasDeuterium<TSpeciesEnum>::value;

template <typename TSpeciesEnum>
inline constexpr bool hasTritium = detail::HasTritium<TSpeciesEnum>::value;
} // namespace psi

template <typename TSpeciesEnum>
class PSIReactionNetwork;

template <typename TSpeciesEnum>
class PSIProductionReaction;

template <typename TSpeciesEnum>
class PSIDissociationReaction;

template <typename TSpeciesEnum>
class PSISinkReaction;
template <typename TSpeciesEnum>
class PSITrapMutationReaction;
template <typename TSpeciesEnum>
class PSIClusterGenerator;

template <typename>
struct TrapMutationClusterData;

enum class PSIFullSpeciesList
{
	He,
	D,
	T,
	V,
	I
};

namespace psi
{
namespace detail
{
template <>
struct HasDeuterium<PSIFullSpeciesList> : std::true_type
{
};

template <>
struct HasTritium<PSIFullSpeciesList> : std::true_type
{
};
} // namespace detail
} // namespace psi

inline const std::string&
toLabelString(PSIFullSpeciesList species)
{
	static const std::string labelArray[] = {"He", "D", "T", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(PSIFullSpeciesList species)
{
	static const std::string nameArray[] = {
		"Helium", "Deuterium", "Tritium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<PSIFullSpeciesList> :
	std::integral_constant<std::size_t, 5>
{
};

template <>
struct NumberOfInterstitialSpecies<PSIFullSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<PSIFullSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct SpeciesForGrouping<PSIFullSpeciesList, 5>
{
	using Sequence = EnumSequence<PSIFullSpeciesList, 5>;
	static constexpr auto first = Sequence(PSIFullSpeciesList::He);
	static constexpr auto last = Sequence(PSIFullSpeciesList::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<PSIFullSpeciesList>
	mapToMomentId(EnumSequence<PSIFullSpeciesList, 5> value)
	{
		if (value == PSIFullSpeciesList::I)
			return 3;
		return value();
	}
};

enum class PSIDeuteriumSpeciesList
{
	He,
	D,
	V,
	I
};

namespace psi
{
namespace detail
{
template <>
struct HasDeuterium<PSIDeuteriumSpeciesList> : std::true_type
{
};
} // namespace detail
} // namespace psi

inline const std::string&
toLabelString(PSIDeuteriumSpeciesList species)
{
	static const std::string labelArray[] = {"He", "D", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(PSIDeuteriumSpeciesList species)
{
	static const std::string nameArray[] = {
		"Helium", "Deuterium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<PSIDeuteriumSpeciesList> :
	std::integral_constant<std::size_t, 4>
{
};

template <>
struct NumberOfInterstitialSpecies<PSIDeuteriumSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<PSIDeuteriumSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct SpeciesForGrouping<PSIDeuteriumSpeciesList, 4>
{
	using Sequence = EnumSequence<PSIDeuteriumSpeciesList, 4>;
	static constexpr auto first = Sequence(PSIDeuteriumSpeciesList::He);
	static constexpr auto last = Sequence(PSIDeuteriumSpeciesList::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<PSIDeuteriumSpeciesList>
	mapToMomentId(EnumSequence<PSIDeuteriumSpeciesList, 4> value)
	{
		if (value == PSIDeuteriumSpeciesList::I)
			return 2;
		return value();
	}
};

enum class PSITritiumSpeciesList
{
	He,
	T,
	V,
	I
};

namespace psi
{
namespace detail
{
template <>
struct HasTritium<PSITritiumSpeciesList> : std::true_type
{
};
} // namespace detail
} // namespace psi

inline const std::string&
toLabelString(PSITritiumSpeciesList species)
{
	static const std::string labelArray[] = {"He", "T", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(PSITritiumSpeciesList species)
{
	static const std::string nameArray[] = {
		"Helium", "Tritium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<PSITritiumSpeciesList> :
	std::integral_constant<std::size_t, 4>
{
};

template <>
struct NumberOfInterstitialSpecies<PSITritiumSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<PSITritiumSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct SpeciesForGrouping<PSITritiumSpeciesList, 4>
{
	using Sequence = EnumSequence<PSITritiumSpeciesList, 4>;
	static constexpr auto first = Sequence(PSITritiumSpeciesList::He);
	static constexpr auto last = Sequence(PSITritiumSpeciesList::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<PSITritiumSpeciesList>
	mapToMomentId(EnumSequence<PSITritiumSpeciesList, 4> value)
	{
		if (value == PSITritiumSpeciesList::I)
			return 2;
		return value();
	}
};

enum class PSIHeliumSpeciesList
{
	He,
	V,
	I
};

inline const std::string&
toLabelString(PSIHeliumSpeciesList species)
{
	static const std::string labelArray[] = {"He", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(PSIHeliumSpeciesList species)
{
	static const std::string nameArray[] = {
		"Helium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<PSIHeliumSpeciesList> :
	std::integral_constant<std::size_t, 3>
{
};

template <>
struct NumberOfInterstitialSpecies<PSIHeliumSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<PSIHeliumSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct SpeciesForGrouping<PSIHeliumSpeciesList, 3>
{
	using Sequence = EnumSequence<PSIHeliumSpeciesList, 3>;
	static constexpr auto first = Sequence(PSIHeliumSpeciesList::He);
	static constexpr auto last = Sequence(PSIHeliumSpeciesList::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<PSIHeliumSpeciesList>
	mapToMomentId(EnumSequence<PSIHeliumSpeciesList, 3> value)
	{
		if (value == PSIHeliumSpeciesList::I)
			return 1;
		return value();
	}
};

template <typename TSpeciesEnum>
struct ReactionNetworkTraits<PSIReactionNetwork<TSpeciesEnum>>
{
	using Species = TSpeciesEnum;

	static constexpr std::size_t numSpecies = numberOfSpecies<TSpeciesEnum>();

	using ProductionReactionType = PSIProductionReaction<Species>;
	using DissociationReactionType = PSIDissociationReaction<Species>;
	using SinkReactionType = PSISinkReaction<Species>;
	using TrapMutationReactionType = PSITrapMutationReaction<Species>;

	using ReactionTypeList = std::tuple<ProductionReactionType,
		DissociationReactionType, SinkReactionType, TrapMutationReactionType>;

	using ClusterGenerator = PSIClusterGenerator<Species>;
};

namespace detail
{
template <typename TSpeciesEnum, typename MemSpace>
struct ClusterDataExtra<PSIReactionNetwork<TSpeciesEnum>, MemSpace>
{
	using NetworkType = PSIReactionNetwork<TSpeciesEnum>;

	ClusterDataExtra() = default;

	template <typename MS>
	KOKKOS_INLINE_FUNCTION
	ClusterDataExtra(const ClusterDataExtra<NetworkType, MS>& data) :
		trapMutationData(data.trapMutationData)
	{
	}

	template <typename MS>
	void
	deepCopy(const ClusterDataExtra<NetworkType, MS>& data)
	{
		trapMutationData.deepCopy(data.trapMutationData);
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return trapMutationData.getDeviceMemorySize();
	}

	using TrapMutationData =
		TrapMutationClusterData<ClusterDataCommon<MemSpace>>;
	TrapMutationData trapMutationData;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
