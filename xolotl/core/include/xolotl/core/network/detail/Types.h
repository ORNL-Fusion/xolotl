#pragma once

#include <xolotl/config.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
using ReactionNetworkIndexType = ::xolotl::IdType;

using CompositionAmountType = ::xolotl::AmountType;

inline constexpr auto invalidNetworkIndex =
	plsm::invalid<ReactionNetworkIndexType>;

inline constexpr auto invalidSpeciesAmount =
	plsm::invalid<CompositionAmountType>;
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
