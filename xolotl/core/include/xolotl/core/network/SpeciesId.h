#pragma once

#include <stdexcept>
#include <type_traits>

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

	operator bool() const noexcept
	{
		return (_id >= 0) && (_id < _nSpec);
	}

	int
	operator()() const noexcept
	{
		return _id;
	}

	template <typename TSpeciesEnum>
	TSpeciesEnum
	cast()
	{
		if (!(*this)) {
			throw InvalidSpeciesId();
		}
		return static_cast<TSpeciesEnum>(_id);
	}

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
