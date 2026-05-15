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
class PSITrapReaction;

template <typename TSpeciesEnum>
class PSITrapMutationReaction;

template <typename TSpeciesEnum>
class PSIClusterGenerator;

template <typename>
struct TrapMutationClusterData;

enum class PSIFullSpeciesList
{
	Trap,
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
	static const std::string labelArray[] = {"Trap", "He", "D", "T", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(PSIFullSpeciesList species)
{
	static const std::string nameArray[] = {
		"Trap", "Helium", "Deuterium", "Tritium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<PSIFullSpeciesList> :
	std::integral_constant<std::size_t, 6>
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
struct SpeciesForGrouping<PSIFullSpeciesList, 6>
{
	using Sequence = EnumSequence<PSIFullSpeciesList, 6>;
	static constexpr auto first = Sequence(PSIFullSpeciesList::He);
	static constexpr auto last = Sequence(PSIFullSpeciesList::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<PSIFullSpeciesList>
	mapToMomentId(EnumSequence<PSIFullSpeciesList, 6> value)
	{
		if (value == PSIFullSpeciesList::I)
			return 3;
		return value() - 1;
	}
};

enum class PSIDeuteriumSpeciesList
{
	Trap,
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
	static const std::string labelArray[] = {"Trap", "He", "D", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(PSIDeuteriumSpeciesList species)
{
	static const std::string nameArray[] = {
		"Trap", "Helium", "Deuterium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<PSIDeuteriumSpeciesList> :
	std::integral_constant<std::size_t, 5>
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
struct SpeciesForGrouping<PSIDeuteriumSpeciesList, 5>
{
	using Sequence = EnumSequence<PSIDeuteriumSpeciesList, 5>;
	static constexpr auto first = Sequence(PSIDeuteriumSpeciesList::He);
	static constexpr auto last = Sequence(PSIDeuteriumSpeciesList::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<PSIDeuteriumSpeciesList>
	mapToMomentId(EnumSequence<PSIDeuteriumSpeciesList, 5> value)
	{
		if (value == PSIDeuteriumSpeciesList::I)
			return 2;
		return value() - 1;
	}
};

enum class PSITritiumSpeciesList
{
	Trap,
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
	static const std::string labelArray[] = {"Trap", "He", "T", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(PSITritiumSpeciesList species)
{
	static const std::string nameArray[] = {
		"Trap", "Helium", "Tritium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<PSITritiumSpeciesList> :
	std::integral_constant<std::size_t, 5>
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
struct SpeciesForGrouping<PSITritiumSpeciesList, 5>
{
	using Sequence = EnumSequence<PSITritiumSpeciesList, 5>;
	static constexpr auto first = Sequence(PSITritiumSpeciesList::He);
	static constexpr auto last = Sequence(PSITritiumSpeciesList::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<PSITritiumSpeciesList>
	mapToMomentId(EnumSequence<PSITritiumSpeciesList, 5> value)
	{
		if (value == PSITritiumSpeciesList::I)
			return 2;
		return value() - 1;
	}
};

enum class PSIHeliumSpeciesList
{
	Trap,
	He,
	V,
	I
};

inline const std::string&
toLabelString(PSIHeliumSpeciesList species)
{
	static const std::string labelArray[] = {"Trap", "He", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(PSIHeliumSpeciesList species)
{
	static const std::string nameArray[] = {
		"Trap", "Helium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<PSIHeliumSpeciesList> :
	std::integral_constant<std::size_t, 4>
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
struct SpeciesForGrouping<PSIHeliumSpeciesList, 4>
{
	using Sequence = EnumSequence<PSIHeliumSpeciesList, 4>;
	static constexpr auto first = Sequence(PSIHeliumSpeciesList::He);
	static constexpr auto last = Sequence(PSIHeliumSpeciesList::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<PSIHeliumSpeciesList>
	mapToMomentId(EnumSequence<PSIHeliumSpeciesList, 4> value)
	{
		if (value == PSIHeliumSpeciesList::I)
			return 1;
		return value() - 1;
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
	using TrapReactionType = PSITrapReaction<Species>;
	using TrapMutationReactionType = PSITrapMutationReaction<Species>;

	using ReactionTypeList =
		std::tuple<ProductionReactionType, DissociationReactionType,
			SinkReactionType, TrapMutationReactionType, TrapReactionType>;

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
