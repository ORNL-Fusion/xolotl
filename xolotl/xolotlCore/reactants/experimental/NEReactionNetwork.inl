#pragma once

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
KOKKOS_INLINE_FUNCTION
void
NEReactionValidator::operator()(std::size_t i, std::size_t j,
    const Subpaving& subpaving, const UpperTriangle<ClusterSet>& prodSet,
    const UpperTriangle<ClusterSet>& dissSet) const
{
}
}

inline
detail::NEReactionValidator
NEReactionNetwork::getReactionValidator() const noexcept
{
    return {};
}
}
}
