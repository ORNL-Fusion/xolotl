#pragma once

#include <stdexcept>
#include <type_traits>

#include <Kokkos_Macros.hpp>

namespace xolotl
{
namespace core
{
namespace network
{
class InvalidSpeciesId : public std::runtime_error
{
public:
	InvalidSpeciesId() : std::runtime_error("Invalid Species Id")
	{
	}

	using std::runtime_error::runtime_error;
};

/**
 * @brief Class representing species id, where a species is a type of
 * atom and is a dimension in the phase space.
 */
class SpeciesId
{
public:
	SpeciesId() = default;

	SpeciesId(int id, int numberOfSpecies) : _id(id), _nSpec(numberOfSpecies)
	{
	}

	template <typename TSpeciesEnum,
		std::enable_if_t<std::is_enum<TSpeciesEnum>::value, int> = 0>
	SpeciesId(TSpeciesEnum species, int numberOfSpecies) :
		_id(static_cast<std::underlying_type_t<TSpeciesEnum>>(species)),
		_nSpec(numberOfSpecies)
	{
	}

	explicit SpeciesId(int numberOfSpecies) : _nSpec(numberOfSpecies)
	{
	}

    KOKKOS_INLINE_FUNCTION
	operator bool() const noexcept
	{
		return (_id >= 0) && (_id < _nSpec);
	}

    KOKKOS_INLINE_FUNCTION
	int
	operator()() const noexcept
	{
		return _id;
	}

	template <typename TSpeciesEnum>
	TSpeciesEnum
	cast() const
	{
		if (!(*this)) {
			throw InvalidSpeciesId();
		}
		return static_cast<TSpeciesEnum>(_id);
	}

	template <typename TSpeciesEnum>
	KOKKOS_INLINE_FUNCTION
	TSpeciesEnum
	unsafeCast() const
	{
		return static_cast<TSpeciesEnum>(_id);
	}

    KOKKOS_INLINE_FUNCTION
	SpeciesId&
	operator++()
	{
		++_id;
		return *this;
	}

private:
	int _id{0};
	int _nSpec;
};
} // namespace network
} // namespace core
} // namespace xolotl
