#pragma once

#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/detail/TrapMutationClusterData.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class TrapMutationHandler
{
public:
	using AmountType = CompositionAmountType;

	TrapMutationHandler() = default;
	virtual ~TrapMutationHandler() = default;

	virtual std::array<std::vector<AmountType>, 7>
	getAllVSizes() const = 0;

	virtual void
	updateData(double temp) = 0;

	const DesorptionInitializer&
	getDesorptionInitializer() const noexcept
	{
		return _desorp;
	}

	const std::array<double, 7>&
	getDepths() const noexcept
	{
		return _depths;
	}

	const std::array<AmountType, 7>&
	getVacancySizes() const noexcept
	{
		return _vSizes;
	}

protected:
	DesorptionInitializer _desorp{};
	std::array<double, 7> _depths{};
	std::array<AmountType, 7> _vSizes{};
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
