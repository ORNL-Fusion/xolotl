#pragma once

#include <xolotl/core/network/IReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
class IPSIReactionNetwork : public IReactionNetwork
{
public:
	using IReactionNetwork::IReactionNetwork;

	virtual SpeciesId
	getHeliumSpeciesId() const = 0;

	virtual SpeciesId
	getVacancySpeciesId() const = 0;

	virtual SpeciesId
	getInterstitialSpeciesId() const = 0;

	virtual bool
	hasDeuterium() const noexcept = 0;

	virtual bool
	hasTritium() const noexcept = 0;

	virtual double
	getTotalTrappedHeliumConcentration(
		ConcentrationsView concentrations, AmountType minSize) = 0;

	virtual void
	updateBurstingConcs(double* gridPointSolution, double factor,
		std::vector<double>& nBurst) = 0;
};
} // namespace network
} // namespace core
} // namespace xolotl
