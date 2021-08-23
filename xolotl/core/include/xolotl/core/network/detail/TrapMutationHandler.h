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
/**
 * @brief Class implementing the parameters for the near surface trap mutation.
 *
 * _desorp : set of helium size and desorption factor
 * _depths: vector indicating at which depth each reaction is occurring, where
 * the index is the helium size - 1 _vSizes : vector indicating how many
 * vacancies are added to the He cluster by the reaction, where the index is the
 * helium size - 1
 */
class TrapMutationHandler
{
public:
	using AmountType = CompositionAmountType;

	TrapMutationHandler() = default;
	virtual ~TrapMutationHandler() = default;

	/**
	 * @brief Method returning all the possible V sizes for each He cluster
	 * size, independent of temperature. This is needed to correctly initialize
	 * the connectivities and be able to switch reaction sets if the temperature
	 * evolves.
	 */
	virtual std::array<std::vector<AmountType>, 7>
	getAllVSizes() const = 0;

	/**
	 * @brief Method updating all the protected class members
	 * depending on the running conditions (material orientation and
	 * temperature).
	 */
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
