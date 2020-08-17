#pragma once

#include <type_traits>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TEnum, std::size_t N>
struct EnumSequence
{
	using Enum = TEnum;
	static_assert(std::is_enum<Enum>::value, "");
	using Rep = std::underlying_type_t<Enum>;
	// static_assert(std::is_signed<Rep>::value, ""); // maybe ?

	KOKKOS_INLINE_FUNCTION
	constexpr EnumSequence(Enum val) : value{val}
	{
	}

	KOKKOS_INLINE_FUNCTION
	constexpr operator Enum() const noexcept
	{
		return value;
	}

	KOKKOS_INLINE_FUNCTION
	constexpr Rep
	operator()() const noexcept
	{
		return static_cast<Rep>(value);
	}

	KOKKOS_INLINE_FUNCTION
	static constexpr std::size_t
	size() noexcept
	{
		return N;
	}

	KOKKOS_INLINE_FUNCTION
	static constexpr EnumSequence
	first() noexcept
	{
		return static_cast<Enum>(0);
	}

	KOKKOS_INLINE_FUNCTION
	static constexpr EnumSequence
	last() noexcept
	{
		return static_cast<Enum>(N - 1);
	}

	KOKKOS_INLINE_FUNCTION
	static constexpr Enum
	invalid() noexcept
	{
		return static_cast<Enum>(-1);
	}

	// TODO: Should these wrap?
	KOKKOS_INLINE_FUNCTION
	EnumSequence&
	operator++()
	{
		value = static_cast<Enum>(static_cast<Rep>(value) + 1);
		return *this;
	}

	KOKKOS_INLINE_FUNCTION
	EnumSequence
	operator++(int)
	{
        EnumSequence ret = *this;
		++static_cast<Rep&>(value);
		return ret;
	}

	Enum value;
};

template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr EnumSequence<T, N>
operator+(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	using Rep = typename EnumSequence<T, N>::Rep;
	return static_cast<Rep>(a.value) + static_cast<Rep>(b.value);
}

template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr EnumSequence<T, N>
operator+(EnumSequence<T, N> a, typename EnumSequence<T, N>::Rep b)
{
	using Rep = typename EnumSequence<T, N>::Rep;
	return static_cast<T>(static_cast<Rep>(a.value) + b);
}

template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr EnumSequence<T, N>
operator-(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	using Rep = typename EnumSequence<T, N>::Rep;
	return static_cast<T>(
		static_cast<Rep>(a.value) - static_cast<Rep>(b.value));
}

template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr EnumSequence<T, N>
operator-(EnumSequence<T, N> a, typename EnumSequence<T, N>::Rep b)
{
	using Rep = typename EnumSequence<T, N>::Rep;
	return static_cast<T>(static_cast<Rep>(a.value) - b);
}

template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr bool
operator<(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	return a.value < b.value;
}

template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr bool
operator<=(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	return a.value <= b.value;
}

template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr bool
operator>(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	return a.value > b.value;
}

template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr bool
operator>=(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	return a.value >= b.value;
}

template <typename TEnum, std::size_t N>
class EnumSequenceRange
{
public:
	using Enum = TEnum;
	using Sequence = EnumSequence<Enum, N>;

	class Iterator
	{
	public:
		// using difference_type = Sequence;
		using value_type = Sequence;
		using reference = const Sequence&;

		KOKKOS_INLINE_FUNCTION
		explicit Iterator(Sequence seq) noexcept : _curr{seq}
		{
		}

		KOKKOS_INLINE_FUNCTION
		constexpr reference
		operator*() const noexcept
		{
			return _curr;
		}

		KOKKOS_INLINE_FUNCTION
		Iterator&
		operator++() noexcept
		{
			++_curr;
			return *this;
		}

		KOKKOS_INLINE_FUNCTION
		constexpr bool
		operator==(Iterator other) const noexcept
		{
			return _curr == other._curr;
		}

		KOKKOS_INLINE_FUNCTION
		constexpr bool
		operator!=(Iterator other) const noexcept
		{
			return !((*this) == other);
		}

	private:
		Sequence _curr;
	};

	KOKKOS_INLINE_FUNCTION
	constexpr EnumSequenceRange()
	{
	}

	KOKKOS_INLINE_FUNCTION
	explicit constexpr EnumSequenceRange(Sequence first, Sequence last) :
		_first{first}, _last{last}
	{
	}

	KOKKOS_INLINE_FUNCTION
	constexpr Iterator
	begin() const noexcept
	{
		return Iterator{_first};
	}

	KOKKOS_INLINE_FUNCTION
	constexpr Iterator
	end() const noexcept
	{
		return Iterator{_last + 1};
	}

private:
	Sequence _first{Sequence::first()};
	Sequence _last{Sequence::last()};
};
} // namespace network
} // namespace core
} // namespace xolotl
