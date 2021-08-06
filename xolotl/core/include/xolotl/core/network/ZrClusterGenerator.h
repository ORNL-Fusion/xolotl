#pragma once

#include <plsm/refine/Detector.h>

namespace xolotl
{
namespace core
{
namespace network
{
class ZrClusterGenerator : public plsm::refine::Detector<ZrClusterGenerator>
{
public:
	using Superclass = plsm::refine::Detector<ZrClusterGenerator>;
	using NetworkType = ZrReactionNetwork;
	using Species = ZrSpecies;

	template <typename PlsmContext>
	using Cluster = Cluster<NetworkType, PlsmContext>;

	using Region = typename NetworkType::Region;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;
	using BoolArray = typename Superclass::BoolVec<Region>;

	ZrClusterGenerator(const options::IOptions& options);

	ZrClusterGenerator(
		const options::IOptions& options, std::size_t refineDepth);

	KOKKOS_INLINE_FUNCTION
	bool
	refine(const Region& region, BoolArray& result) const;

	KOKKOS_INLINE_FUNCTION
	bool
	select(const Region& region) const;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getFormationEnergy(const Cluster<PlsmContext>& cluster) const noexcept;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getMigrationEnergy(const Cluster<PlsmContext>& cluster) const noexcept;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getDiffusionFactor(const Cluster<PlsmContext>& cluster,
		double latticeParameter) const noexcept;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getReactionRadius(const Cluster<PlsmContext>& cluster,
		double latticeParameter, double interstitialBias,
		double impurityRadius) const noexcept;

private:
	AmountType _maxV;
	AmountType _maxI;
};
} // namespace network
} // namespace core
} // namespace xolotl
