#pragma once

#include <type_traits>

namespace xolotl
{
namespace core
{
namespace network
{
/**
 * @brief Treat the set of values of an enum as a sequence
 *
 * This class wraps a single enum value and provides operations as if the
 * set of enum values forms a sequence.
 *
 * @note: This class assumes that the enum values have the default
 * representation numbering. That is, it is assumed that the underlying
 * integral type values form a sequence of values beginning at zero.
 *
 * @tparam TEnum The enum type (scoped or unscoped)
 * @tparam N The number of valid values in the enum
 */
template <typename TEnum, std::size_t N>
struct EnumSequence
{
	static_assert(std::is_enum_v<TEnum>, "");
	using Enum = TEnum;
	using Rep = std::underlying_type_t<Enum>;
	// static_assert(std::is_signed<Rep>::value, ""); // maybe ?

	/**
	 * @brief Convert from enum type
	 */
	KOKKOS_INLINE_FUNCTION
	constexpr EnumSequence(Enum val) : value{val}
	{
	}

	/**
	 * @brief Convert to enum type
	 */
	KOKKOS_INLINE_FUNCTION
	constexpr operator Enum() const noexcept
	{
		return value;
	}

	/**
	 * @brief Get integral representation value
	 */
	KOKKOS_INLINE_FUNCTION
	constexpr Rep
	operator()() const noexcept
	{
		return static_cast<Rep>(value);
	}

	/**
	 * @brief Get number of enum values
	 */
	KOKKOS_INLINE_FUNCTION
	static constexpr std::size_t
	size() noexcept
	{
		return N;
	}

	/**
	 * @brief Get first enum value
	 */
	KOKKOS_INLINE_FUNCTION
	static constexpr EnumSequence
	first() noexcept
	{
		return static_cast<Enum>(0);
	}

	/**
	 * @brief Get last enum value
	 */
	KOKKOS_INLINE_FUNCTION
	static constexpr EnumSequence
	last() noexcept
	{
		return static_cast<Enum>(N - 1);
	}

	/**
	 * @brief Get a value that is not a valid enum value
	 */
	KOKKOS_INLINE_FUNCTION
	static constexpr Enum
	invalid() noexcept
	{
		return static_cast<Enum>(-1);
	}

	/**
	 * @brief Pre-increment
	 * @todo: Should these wrap?
	 */
	KOKKOS_INLINE_FUNCTION
	EnumSequence&
	operator++()
	{
		value = static_cast<Enum>(static_cast<Rep>(value) + 1);
		return *this;
	}

	/**
	 * @brief Post-increment
	 */
	KOKKOS_INLINE_FUNCTION
	EnumSequence
	operator++(int)
	{
		EnumSequence ret = *this;
		++static_cast<Rep&>(value);
		return ret;
	}

	/** Wrapped enum value */
	Enum value;
};

/**
 * @relates EnumSequence
 * @brief Addition operator
 */
template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr EnumSequence<T, N>
operator+(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	using Rep = typename EnumSequence<T, N>::Rep;
	return static_cast<Rep>(a.value) + static_cast<Rep>(b.value);
}

/**
 * @relates EnumSequence
 * @brief Addition operator
 */
template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr EnumSequence<T, N>
operator+(EnumSequence<T, N> a, typename EnumSequence<T, N>::Rep b)
{
	using Rep = typename EnumSequence<T, N>::Rep;
	return static_cast<T>(static_cast<Rep>(a.value) + b);
}

/**
 * @relates EnumSequence
 * @brief Subtraction operator
 */
template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr EnumSequence<T, N>
operator-(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	using Rep = typename EnumSequence<T, N>::Rep;
	return static_cast<T>(
		static_cast<Rep>(a.value) - static_cast<Rep>(b.value));
}

/**
 * @relates EnumSequence
 * @brief Subtraction operator
 */
template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr EnumSequence<T, N>
operator-(EnumSequence<T, N> a, typename EnumSequence<T, N>::Rep b)
{
	using Rep = typename EnumSequence<T, N>::Rep;
	return static_cast<T>(static_cast<Rep>(a.value) - b);
}

/**
 * @relates EnumSequence
 * @brief Comparison operator
 */
template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr bool
operator<(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	return a.value < b.value;
}

/**
 * @relates EnumSequence
 * @brief Comparison operator
 */
template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr bool
operator<=(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	return a.value <= b.value;
}

/**
 * @relates EnumSequence
 * @brief Comparison operator
 */
template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr bool
operator>(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	return a.value > b.value;
}

/**
 * @relates EnumSequence
 * @brief Comparison operator
 */
template <typename T, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr bool
operator>=(EnumSequence<T, N> a, EnumSequence<T, N> b)
{
	return a.value >= b.value;
}

/**
 * @brief Represent an iterable range of values in an enum
 *
 * This class enables iteration over the values of an enum in a
 * range-based for loop
 *
 * @note Makes same assumptions as EnumSequence
 *
 * @tparam TEnum The enum type (scoped or unscoped)
 * @tparam N The number of valid values in the enum
 */
template <typename TEnum, std::size_t N>
class EnumSequenceRange
{
public:
	using Enum = TEnum;
	using Sequence = EnumSequence<Enum, N>;

	/**
	 * @brief Implements standard iterator interface
	 */
	class Iterator
	{
	public:
		// using difference_type = Sequence;
		using value_type = Sequence;
		using reference = const Sequence&;

		/**
		 * @brief Construct from sequence value
		 */
		KOKKOS_INLINE_FUNCTION
		explicit Iterator(Sequence seq) noexcept : _curr{seq}
		{
		}

		/**
		 * @brief Dereference
		 */
		KOKKOS_INLINE_FUNCTION
		constexpr reference
		operator*() const noexcept
		{
			return _curr;
		}

		/**
		 * @brief Increment
		 */
		KOKKOS_INLINE_FUNCTION
		Iterator&
		operator++() noexcept
		{
			++_curr;
			return *this;
		}

		/**
		 * @brief Comparison operator
		 */
		KOKKOS_INLINE_FUNCTION
		constexpr bool
		operator==(Iterator other) const noexcept
		{
			return _curr == other._curr;
		}

		/**
		 * @brief Comparison operator
		 */
		KOKKOS_INLINE_FUNCTION
		constexpr bool
		operator!=(Iterator other) const noexcept
		{
			return !((*this) == other);
		}

	private:
		Sequence _curr;
	};

	/**
	 * @brief Default construct uses first and last sequence values
	 */
	KOKKOS_INLINE_FUNCTION
	constexpr EnumSequenceRange()
	{
	}

	/**
	 * @brief Construct with custom range endpoints
	 */
	KOKKOS_INLINE_FUNCTION
	explicit constexpr EnumSequenceRange(Sequence first, Sequence last) :
		_first{first},
		_last{last}
	{
	}

	/**
	 * @brief Get iterator to first value
	 */
	KOKKOS_INLINE_FUNCTION
	constexpr Iterator
	begin() const noexcept
	{
		return Iterator{_first};
	}

	/**
	 * @brief Get iterator to one-past-the-last value
	 */
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
