#pragma once

#include <xolotl/core/network/ReactionNetwork.h>
#include <xolotl/core/network/T91Reaction.h>
#include <xolotl/core/network/T91Traits.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class T91ReactionGenerator;
}

class T91ReactionNetwork : public ReactionNetwork<T91ReactionNetwork>
{
	friend class ReactionNetwork<T91ReactionNetwork>;

public:
	using Superclass = ReactionNetwork<T91ReactionNetwork>;
	using Subpaving = typename Superclass::Subpaving;
	using Composition = typename Superclass::Composition;
	using Species = typename Superclass::Species;
	using AmountType = typename Superclass::AmountType;
	using IndexType = typename Superclass::IndexType;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;

	using Superclass::Superclass;

	IndexType
	checkLargestClusterId();

	std::string
	getMonitorOutputFileName() const override
	{
		return "T91.dat";
	}

	std::string
	getMonitorDataHeaderString() const override;

	void
	addMonitorDataValues(Kokkos::View<const double*> conc, double fac,
		std::vector<double>& totalVals) override;

	std::size_t
	getMonitorDataLineSize() const override
	{
		return getSpeciesListSize() * 4;
	}

	void
	writeMonitorDataLine(
		const std::vector<double>& localData, double time) override;

private:
	double
	checkLatticeParameter(double latticeParameter);

	double
	computeAtomicVolume(double latticeParameter)
	{
		// 2 atoms per cell
		return 0.5 * latticeParameter * latticeParameter * latticeParameter;
	}

	double
	checkImpurityRadius(double impurityRadius);

	detail::T91ReactionGenerator
	getReactionGenerator() const noexcept;

	void
	readClusters(const std::string filename)
	{
		return;
	}

	void
	readReactions(double temperature, const std::string filename)
	{
		return;
	}
};

namespace detail
{
class T91ReactionGenerator :
	public ReactionGenerator<T91ReactionNetwork, T91ReactionGenerator>
{
	friend class ReactionGeneratorBase<T91ReactionNetwork,
		T91ReactionGenerator>;

public:
	using NetworkType = T91ReactionNetwork;
	using Subpaving = typename NetworkType::Subpaving;
	using IndexType = typename NetworkType::IndexType;

	using Superclass =
		ReactionGenerator<T91ReactionNetwork, T91ReactionGenerator>;

	using Superclass::Superclass;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	operator()(IndexType i, IndexType j, TTag tag) const;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	addSinks(IndexType i, TTag tag) const;

private:
	ReactionCollection<NetworkType>
	getReactionCollection() const;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/T91ClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/T91ReactionNetwork.tpp>
#endif
